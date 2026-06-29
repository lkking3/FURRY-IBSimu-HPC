// ============================================================================
//  meniscus_cell.cpp
//  High-resolution AXISYMMETRIC (MODE_CYL) self-consistent meniscus / extraction
//  micro-model for ONE bore of the curved two-grid system.
//
//  Why: the plasma meniscus sheath is ~Debye-length thick (~20 um at n_i~4e17),
//  far finer than the full 3-D grid mesh (0.2-0.3 mm). It cannot be resolved in
//  the full model, so we solve it here in a small fine-mesh cell and export the
//  extracted beamlet phase space, which the 3-D model injects per aperture
//  (see MENISCUS_UPGRADE_PLAN.md and grid_3d_curved.cpp INJECTION_MODE=meniscus).
//
//  Geometry (axisymmetric, axial = x, radial = r = y >= 0): a plasma region
//  upstream, a screen electrode (V=VS) and an accel electrode (V=VA), each a flat
//  slab pierced by an on-axis bore. Local curvature over one ~3 mm bore is
//  negligible, so a flat-gap cell is an excellent local approximation; the
//  beamlet's aim (tilt toward the focus) is applied later when it is placed in 3-D.
//
//  Output: results_men/meniscus_beamlet.dat  -- (r, r', current) at the accel
//  exit plane, plus a header (energy, total current per bore, n points).
//
//  Build: make -f Makefile.3d meniscus_cell      Run via SLURM (single bore is
//  cheap, but keep all compute off the login node).
// ============================================================================
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <sys/stat.h>

#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_bicgstabsolver.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "meshscalarfield.hpp"
#include "particledatabase.hpp"      // ParticleDataBaseCyl
#include "trajectorydiagnostics.hpp"
#include "config.h"                  // InitialPlasma, AXIS_*
#include "error.hpp"
#include "ibsimu.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double envd(const char* k,double v){const char*s=getenv(k);return s?strtod(s,nullptr):v;}
static int    envi(const char* k,int v){const char*s=getenv(k);return s?atoi(s):v;}
static std::string envs(const char* k,const char* v){const char*s=getenv(k);return s?std::string(s):std::string(v);}
static bool file_ok(const std::string &p){ std::ifstream f(p.c_str()); return f.good(); }
static void mkdirs(const std::string &dir){
    std::string acc;
    for( size_t i=0;i<dir.size();++i ){ if(dir[i]=='/'&&!acc.empty()) mkdir(acc.c_str(),0777); acc+=dir[i]; }
    if( !acc.empty() ) mkdir(acc.c_str(),0777);
}
static void mkdir_parent(const std::string &p){ size_t s=p.find_last_of('/'); if(s!=std::string::npos) mkdirs(p.substr(0,s)); }

static const double QE  = 1.602176634e-19;
static const double AMU = 1.66053906660e-27;

// --- electrode geometry globals (set in simu, read by FuncSolid callbacks) ---
static double g_xs0, g_xs1, g_rbs;   // screen: x in [xs0,xs1], material at r >= rbs
static double g_xa0, g_xa1, g_rba;   // accel : x in [xa0,xa1], material at r >= rba
static bool screen_fn(double x,double r,double /*z*/){ return x>=g_xs0 && x<=g_xs1 && r>=g_rbs; }
static bool accel_fn (double x,double r,double /*z*/){ return x>=g_xa0 && x<=g_xa1 && r>=g_rba; }

void simu( int argc, char **argv )
{
    ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
    ibsimu.set_thread_count( envi("ION_THREADS", 1) );

    // ---- geometry parameters (SI metres) ----
    const double H        = envd("MEN_H", 5.0e-6);          // fine mesh (resolve Debye)
    const double T_SCR    = envd("MEN_T_SCR", 6.0e-3);
    const double T_ACC    = envd("MEN_T_ACC", 6.0e-3);
    const double GAP      = envd("MEN_GAP", 4.0e-3);
    const double RB_SCR   = envd("MEN_RBORE_SCR", 1.5e-3);  // screen bore radius
    const double RB_ACC   = envd("MEN_RBORE_ACC", 1.5e-3);  // accel  bore radius
    const double RMAX     = envd("MEN_RMAX", 4.0e-3);       // radial domain
    const double X_PLASMA = envd("MEN_X_PLASMA", 0.5e-3);   // plasma fills x < this
    const double DRIFT    = envd("MEN_DRIFT", 4.0e-3);      // downstream of accel

    const double VS = envd("VS_V", 0.0);
    const double VA = envd("VA_V", -10000.0);

    // ---- plasma + species (match the 3-D run) ----
    const double NI    = envd("PLASMA_NI_M3", 4.0e17);
    const double TE    = envd("PLASMA_TE_EV", 3.7);
    const double Q_E   = envd("ION_Q_E", 1.0);
    const double M_AMU = envd("ION_M_AMU", 2.0);
    const double TP    = envd("ION_TP_EV", 0.0);
    const double TT    = envd("ION_TT_EV", 0.2);
    const double JSC   = envd("ION_J_SCALE", 1.0);
    const int    NPART = envi("MEN_NPART", 4000);
    const int    ITERM = envi("MEN_ITER_MAX", 18);

    const std::string OUTDIR = envs("MEN_RESULTS_DIR","results_men");
    mkdirs( OUTDIR );
    auto opath = [&](const char* f){ return OUTDIR + "/" + f; };

    // Beamlet output path (in its own cache subfolder) + reuse cache: if a solved
    // beamlet already exists here and MENISCUS_REUSE!=0, skip the (cheap) solve.
    // Delete the file or set MENISCUS_REUSE=0 to re-solve after changing params.
    const std::string OUT_BEAM = envs("MENISCUS_FILE","meniscus_cache/meniscus_beamlet.dat");
    const int         REUSE    = envi("MENISCUS_REUSE", 1);
    if( REUSE && file_ok(OUT_BEAM) ) {
        std::printf("meniscus: reusing cached beamlet '%s' (set MENISCUS_REUSE=0 to force re-solve)\n",
                    OUT_BEAM.c_str());
        return;
    }

    // ---- axial layout: plasma | screen | gap | accel | drift ----
    g_xs0 = X_PLASMA;                 g_xs1 = g_xs0 + T_SCR;   g_rbs = RB_SCR;
    g_xa0 = g_xs1 + GAP;              g_xa1 = g_xa0 + T_ACC;    g_rba = RB_ACC;
    const double LX     = g_xa1 + DRIFT;
    const double x_exit = g_xa1 + 0.25*DRIFT;     // beamlet sampled just past accel
    const double x_emit = std::max( 2.0*H, g_xs0 - 2.0*H );   // emit just upstream of screen
    const double r_emit = std::min( RB_SCR, RMAX - 2.0*H );

    Int3D  size( (int)floor(LX/H)+1, (int)floor(RMAX/H)+1, 1 );
    Geometry geom( MODE_CYL, size, Vec3D(0.0, 0.0, 0.0), H );
    std::printf("meniscus cell: axial %.2f mm x radial %.2f mm, h=%.1f um -> %d x %d nodes\n",
                LX*1e3, RMAX*1e3, H*1e6, size[0], size[1]);

    geom.set_solid( 7, new FuncSolid(&screen_fn) );
    geom.set_solid( 8, new FuncSolid(&accel_fn) );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET, VS) );   // xmin: plasma chamber
    geom.set_boundary( 2, Bound(BOUND_NEUMANN,   0.0) );  // xmax: downstream
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,   0.0) );  // r=0 axis
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,   0.0) );  // r=RMAX outer
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET, VS) );   // screen
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, VA) );   // accel
    geom.build_mesh();

    EpotField       epot( geom );
    MeshScalarField scharge( geom );
    MeshVectorField bfield;
    EpotEfield      efield( epot );
    EpotBiCGSTABSolver solver( geom );

    InitialPlasma initp( AXIS_X, g_xs0 );
    solver.set_initial_plasma( VS, &initp );
    const double rhoe = -QE * NI;
    solver.set_pexp_plasma( rhoe, TE, VS );

    ParticleDataBaseCyl pdb( geom );
    bool pmirror[6] = { false,false,false,false,false,false };
    pdb.set_mirror( pmirror );

    const double Te_J = TE*QE, mi = M_AMU*AMU;
    const double cs   = std::sqrt( std::max(Te_J/mi, 0.0) );
    const double J    = JSC * 0.6 * Q_E * QE * NI * cs;     // Bohm flux [A/m^2]
    const double E0   = std::max( 1.0, std::fabs(VS - VA) );

    std::printf("J=%.4e A/m^2  E0=%.1f eV  species q=%g m=%g amu  bore r_s=%.2f mm\n",
                J, E0, Q_E, M_AMU, RB_SCR*1e3);

    // ---- self-consistent Vlasov iteration (meniscus forms over iterations) ----
    std::vector<double> hist_I;
    for( int it=0; it<ITERM; ++it ) {
        solver.solve( epot, scharge );
        efield.recalculate();
        pdb.clear();
        pdb.add_2d_beam_with_energy( (uint32_t)NPART, J, Q_E, M_AMU, E0, TP, TT,
                                     x_emit, 0.0, x_emit, r_emit );
        pdb.iterate_trajectories( scharge, efield, bfield );

        // extracted current at the exit plane (convergence monitor)
        std::vector<trajectory_diagnostic_e> dg;
        dg.push_back(DIAG_CURR);
        TrajectoryDiagnosticData td;
        double Iexit = 0.0;
        try {
            pdb.trajectories_at_plane( td, AXIS_X, x_exit, dg );
            const std::vector<double>& c = td(0).data();
            for( double v : c ) Iexit += v;
        } catch(...) {}
        hist_I.push_back( Iexit );
        std::printf("[men iter %d/%d] I_exit=%.4e A  npart=%u\n",
                    it+1, ITERM, Iexit, (unsigned)pdb.size());
        std::fflush(stdout);
        if( it>=4 && hist_I.size()>=2 ) {
            double d = std::fabs(hist_I.back()-hist_I[hist_I.size()-2])
                       / std::max(1e-30, std::fabs(hist_I.back()));
            if( d < 0.005 ) { std::printf("[men] converged (dI=%.2f%%)\n",100*d); break; }
        }
    }

    // ---- export beamlet phase space at the accel exit ----
    std::vector<trajectory_diagnostic_e> diag;
    diag.push_back(DIAG_Y);     // radial position r
    diag.push_back(DIAG_YP);    // radial slope r' = dr/dx
    diag.push_back(DIAG_CURR);  // current carried by trajectory [A]
    TrajectoryDiagnosticData tdata;
    pdb.trajectories_at_plane( tdata, AXIS_X, x_exit, diag );
    const std::vector<double>& r  = tdata(0).data();
    const std::vector<double>& rp = tdata(1).data();
    const std::vector<double>& cu = tdata(2).data();
    const size_t N = cu.size();

    double Itot=0.0, sr2=0.0, srp2=0.0, srrp=0.0, rmax=0.0;
    for( size_t i=0;i<N;++i ){
        Itot += cu[i]; sr2 += cu[i]*r[i]*r[i]; srp2 += cu[i]*rp[i]*rp[i];
        srrp += cu[i]*r[i]*rp[i]; if(r[i]>rmax) rmax=r[i];
    }
    double rrms  = Itot>0? std::sqrt(sr2/Itot):0.0;
    double rprms = Itot>0? std::sqrt(srp2/Itot):0.0;
    double corr  = Itot>0? srrp/Itot:0.0;
    double emit  = std::sqrt( std::max(0.0, (sr2/std::max(Itot,1e-30))*(srp2/std::max(Itot,1e-30)) - corr*corr) );

    // ---- current accounting: analytic Bohm vs injected vs extracted ----
    // Resolves whether I_per_bore is physical or a cyl current-normalization
    // artifact. If extracted ~= injected but both ~= 3.x x analytic Bohm, the
    // factor is in how add_2d_beam's J maps to total revolved current (a
    // convention to divide out). If injected ~= analytic Bohm but extracted is
    // larger, that's nonphysical (a dump bug).
    double I_inj = 0.0;
    {   std::vector<trajectory_diagnostic_e> dgi; dgi.push_back(DIAG_CURR);
        TrajectoryDiagnosticData ti;
        try { pdb.trajectories_at_plane( ti, AXIS_X, x_emit + 2.0*H, dgi );
              for( double v : ti(0).data() ) I_inj += v; } catch(...) {}
    }
    const double I_bohm = J * M_PI * r_emit * r_emit;      // analytic Bohm over emission disk
    std::printf("CURRENT CHECK: analytic Bohm(J*pi*r_emit^2)=%.4e A | injected(@emit)=%.4e A | extracted(@exit)=%.4e A\n",
                I_bohm, I_inj, Itot );
    std::printf("  ratios: extracted/injected=%.3f  extracted/Bohm=%.3f  injected/Bohm=%.3f\n",
                (I_inj>0?Itot/I_inj:0.0), (I_bohm>0?Itot/I_bohm:0.0), (I_bohm>0?I_inj/I_bohm:0.0) );

    mkdir_parent( OUT_BEAM );
    std::ofstream f( OUT_BEAM.c_str() );
    f << "# meniscus beamlet (axisymmetric) at accel exit; SI units\n";
    f << "# E_eV " << E0 << "  I_per_bore_A " << Itot << "  n_points " << N << "\n";
    f << "# r_rms_m " << rrms << "  rp_rms_rad " << rprms
      << "  rms_emittance_m_rad " << emit << "  r_max_m " << rmax << "\n";
    f << "# analytic_Bohm_A " << I_bohm << "  injected_A " << I_inj
      << "  extracted_over_Bohm " << (I_bohm>0?Itot/I_bohm:0.0) << "\n";
    f << "# r_m  rp_rad  curr_A\n";
    f.precision(8);
    for( size_t i=0;i<N;++i )
        f << r[i] << " " << rp[i] << " " << cu[i] << "\n";
    f.close();

    geom.save( opath("men_geom.dat") );
    epot.save( opath("men_epot.dat") );
    std::printf("MENISCUS DONE: I_per_bore=%.4e A  r_rms=%.3f mm  rp_rms=%.2f mrad  "
                "emit=%.3e m.rad  (%zu pts) -> %s\n",
                Itot, rrms*1e3, rprms*1e3, emit, N, OUT_BEAM.c_str());
}

int main( int argc, char **argv )
{
    try { simu( argc, argv ); }
    catch( Error e ) { e.print_error_message( ibsimu.message(0) ); return 1; }
    return 0;
}
