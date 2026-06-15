// ============================================================================
//  grid_3d_curved.cpp
//  3-D IBSimu extraction model for a curved two-grid ion-optics system.
//
//  Geometry  : two spherical-cap electrodes imported as STL solids
//              (screen @ V_s, accel @ V_a), assembled concentric with a fixed
//              face-to-face gap.  The model is a HALF of the grid: the x = 0
//              plane is a Cartesian mirror (Neumann) symmetry plane -- the
//              valid reduction for the hexagonal (6-fold) aperture pattern,
//              which admits NO orthogonal mirror pair (a quarter is invalid).
//
//  Physics   : self-consistent (Vlasov) space charge.  Positive ions are
//              emitted at the Bohm velocity from each aperture along its local
//              surface normal (add_cylindrical_beam_with_energy).  Electrons
//              are treated with the positive-exponential Boltzmann plasma model
//              (set_pexp_plasma), carried over from the validated 2-D module.
//
//  Units     : SI metres internally.  STL files are millimetres -> scaled by a
//              Transformation (1e-3).  The aperture file is ALREADY in metres.
//
//  Build/run : single-threaded-friendly; configure threads via ION_THREADS.
//              All compute on the cluster via SLURM (never the login node).
//
//  NOTE      : IBSimu is not available in the prep sandbox, so this file is
//              written against the IBSimu 1.0.6 API (RADIS STL example +
//              the project's multi_grid_2d_curved.cpp) and must be COMPILED
//              AND LINKED ON THE CLUSTER, where minor signature drift in the
//              diagnostics block (if any) should be reconciled.
// ============================================================================

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <sys/stat.h>

#include "geometry.hpp"
#include "stl_solid.hpp"
#include "stlfile.hpp"
#include "transformation.hpp"
#include "epot_bicgstabsolver.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "meshscalarfield.hpp"
#include "particledatabase.hpp"      // ParticleDataBase3D
#include "trajectorydiagnostics.hpp"
#include "geomplotter.hpp"
#include "config.h"                  // InitialPlasma, AXIS_*
#include "error.hpp"
#include "ibsimu.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ----------------------------- env helpers ---------------------------------
static double envd(const char* k, double v){ const char* s=getenv(k); return s?strtod(s,nullptr):v; }
static int    envi(const char* k, int v){ const char* s=getenv(k); return s?atoi(s):v; }
static std::string envs(const char* k, const char* v){ const char* s=getenv(k); return s?std::string(s):std::string(v); }

// physical constants
static const double QE  = 1.602176634e-19;   // C
static const double AMU = 1.66053906660e-27; // kg

struct Aperture { double cx,cy,cz, nx,ny,nz, r; };

// Read aperture emission sites (SI metres, sim coords) produced by
// tools/extract_apertures_3d.py
static std::vector<Aperture> read_apertures( const std::string &path )
{
    std::vector<Aperture> out;
    std::ifstream in(path.c_str());
    if( !in.good() )
        throw( Error( ERROR_LOCATION, (std::string)"cannot open aperture file: "+path ) );
    std::string line;
    while( std::getline(in,line) ) {
        if( line.empty() || line[0]=='#' ) continue;
        std::istringstream ss(line);
        Aperture a;
        if( ss >> a.cx >> a.cy >> a.cz >> a.nx >> a.ny >> a.nz >> a.r )
            out.push_back(a);
    }
    return out;
}

// Build an orthonormal pair (d1,d2) spanning the plane perpendicular to n, such
// that d1 x d2 = n (so add_cylindrical_beam_with_energy launches along +n).
static void tangent_basis( const Vec3D &n, Vec3D &d1, Vec3D &d2 )
{
    Vec3D a = ( std::fabs(n[0]) < 0.9 ) ? Vec3D(1,0,0) : Vec3D(0,1,0);
    d1 = cross( a, n ); d1 /= d1.norm2();
    d2 = cross( n, d1 ); d2 /= d2.norm2();   // d1 x d2 = n
}

void simu( int argc, char **argv )
{
    ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
    ibsimu.set_thread_count( envi("ION_THREADS", 1) );

    // ----------------------------- configuration ---------------------------
    const double H      = envd("H", 2.0e-4);            // cell size [m]
    const double LX     = envd("LX", 0.095);            // half radial extent (x>=0) [m]
    const double LY     = envd("LY", 0.095);            // +/- transverse extent [m]
    const double LZ     = envd("LZ", 0.2134);           // axial length [m]

    const std::string SCREEN_STL = envs("SCREEN_STL","geometry/screen_half_pos.stl");
    const std::string ACCEL_STL  = envs("ACCEL_STL", "geometry/accel_half_pos.stl");
    const std::string AP_FILE    = envs("APERTURE_FILE","geometry/apertures_screen.dat");
    const double STL_SCALE = envd("STL_SCALE", 1.0e-3); // mm -> m

    const double VS = envd("VS_V", 0.0);                // screen / plasma potential
    const double VA = envd("VA_V", -10000.0);           // accel potential
    const double Z_PLASMA = envd("Z_PLASMA", 0.007);    // plasma fills z < Z_PLASMA [m]

    // plasma + ion species (carried over from the 2-D curved module)
    const double NI   = envd("PLASMA_NI_M3", 1.0e16);
    const double TE   = envd("PLASMA_TE_EV", 3.7);
    const double Q_E  = envd("ION_Q_E", 1.0);
    const double M_AMU= envd("ION_M_AMU", 40.0);        // Ar+ default
    const int    NPART= envi("ION_NPART", 60000);
    const int    ITERM= envi("ION_ITER_MAX", 12);
    const double TP   = envd("ION_TP_EV", 0.0);
    const double TT   = envd("ION_TT_EV", 0.2);
    const double JSC  = envd("ION_J_SCALE", 1.0);
    const double SC_ALPHA = envd("SC_ALPHA", 0.5);      // space-charge averaging

    const std::string OUTDIR = envs("RESULTS_DIR","results_3d");
    mkdir( OUTDIR.c_str(), 0777 );
    auto opath = [&](const char* f){ return OUTDIR + "/" + f; };

    // ----------------------------- geometry --------------------------------
    Int3D meshsize( (int)floor(LX/H)+1, (int)floor(2.0*LY/H)+1, (int)floor(LZ/H)+1 );
    Vec3D origo( 0.0, -LY, 0.0 );
    Geometry geom( MODE_3D, meshsize, origo, H );

    std::printf("3D mesh: %d x %d x %d = %.0f M nodes (h=%.3g m)\n",
                meshsize[0], meshsize[1], meshsize[2],
                (double)meshsize[0]*meshsize[1]*meshsize[2]/1e6, H );

    Transformation T;                       // STL (mm) -> sim (m)
    T.scale( Vec3D( STL_SCALE, STL_SCALE, STL_SCALE ) );

    STLFile *fscr = new STLFile( SCREEN_STL );
    STLSolid *screen = new STLSolid;
    screen->set_transformation( T );
    screen->add_stl_file( fscr );
    geom.set_solid( 7, screen );

    STLFile *facc = new STLFile( ACCEL_STL );
    STLSolid *accel = new STLSolid;
    accel->set_transformation( T );
    accel->add_stl_file( facc );
    geom.set_solid( 8, accel );

    // boundary ids: 1 xmin, 2 xmax, 3 ymin, 4 ymax, 5 zmin, 6 zmax
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,   0.0) );  // x=0 mirror (symmetry)
    geom.set_boundary( 2, Bound(BOUND_NEUMANN,   0.0) );  // outer radial
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,   0.0) );  // outer transverse
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,   0.0) );  // outer transverse
    geom.set_boundary( 5, Bound(BOUND_DIRICHLET, VS ) );  // upstream plasma chamber
    geom.set_boundary( 6, Bound(BOUND_NEUMANN,   0.0) );  // downstream drift (open)
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET, VS ) );  // screen
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, VA ) );  // accel

    geom.build_mesh();
    geom.build_surface();

    // ----------------------------- solver / plasma -------------------------
    EpotBiCGSTABSolver solver( geom );
    InitialPlasma initp( AXIS_Z, Z_PLASMA );
    solver.set_initial_plasma( VS, &initp );
    const double rhoe = -QE * NI;                       // quasi-neutral n_e ~ n_i
    solver.set_pexp_plasma( rhoe, TE, VS );

    EpotField       epot( geom );
    MeshScalarField scharge( geom );
    MeshScalarField scharge_ave( geom );
    MeshVectorField bfield;                             // zero B
    EpotEfield      efield( epot );
    field_extrpl_e ex[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                             FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                             FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( ex );

    ParticleDataBase3D pdb( geom );
    pdb.set_max_steps( envi("ION_MAX_STEPS", 4000) );
    bool pmirror[6] = { true, false, false, false, false, false }; // mirror x=0
    pdb.set_mirror( pmirror );

    // ----------------------------- injection setup -------------------------
    std::vector<Aperture> aps = read_apertures( AP_FILE );
    const int n_per = std::max( 1, NPART / std::max((size_t)1, aps.size()) );
    const double Te_J = TE * QE;
    const double mi   = M_AMU * AMU;
    const double cs   = std::sqrt( std::max(Te_J/mi, 0.0) );        // Bohm speed
    const double J    = JSC * 0.6 * Q_E * QE * NI * cs;             // [A/m^2]
    const double E0   = std::max( 1.0, std::fabs(VS - VA) );        // [eV]

    std::printf("apertures=%zu  n_per=%d  J=%.3e A/m^2  E0=%.1f eV  species: q=%g m=%g amu\n",
                aps.size(), n_per, J, E0, Q_E, M_AMU );

    // ----------------------------- Vlasov iteration ------------------------
    for( int it = 0; it < ITERM; ++it ) {
        std::printf("[iter %d/%d] solving Poisson ...\n", it+1, ITERM); std::fflush(stdout);
        solver.solve( epot, scharge_ave );
        efield.recalculate();

        std::printf("[iter %d/%d] injecting %zu apertures ...\n", it+1, ITERM, aps.size());
        pdb.clear();
        for( const Aperture &a : aps ) {
            Vec3D n( a.nx, a.ny, a.nz );
            Vec3D d1, d2; tangent_basis( n, d1, d2 );
            // emit from 2h upstream of the surface, inside the plasma region
            Vec3D c( a.cx - 2.0*H*n[0], a.cy - 2.0*H*n[1], a.cz - 2.0*H*n[2] );
            pdb.add_cylindrical_beam_with_energy(
                (uint32_t)n_per, J, Q_E, M_AMU,
                E0, TP, TT,
                c, d1, d2, a.r );
        }
        pdb.iterate_trajectories( scharge, efield, bfield );
        std::printf("[iter %d/%d] traced npart=%u\n", it+1, ITERM, (unsigned)pdb.size());

        // space-charge averaging for stability (RADIS-style)
        if( it == 0 ) {
            scharge_ave = scharge;
        } else {
            const double b = 1.0 - SC_ALPHA;
            uint32_t nc = scharge.nodecount();
            for( uint32_t k = 0; k < nc; ++k )
                scharge_ave(k) = SC_ALPHA*scharge(k) + b*scharge_ave(k);
        }
    }

    // ----------------------------- save state ------------------------------
    geom.save( opath("geom.dat") );
    epot.save( opath("epot.dat") );
    pdb.save(  opath("pdb.dat")  );

    // ----------------------------- diagnostics -----------------------------
    // Merged-beam divergence at a downstream plane and accel interception.
    // Uses TrajectoryDiagnosticData at a z = const plane.
    try {
        const double z_diag = envd("Z_DIAG", LZ - 5.0*H);
        std::vector<trajectory_diagnostic_e> diag;
        diag.push_back( DIAG_X );  diag.push_back( DIAG_VX );
        diag.push_back( DIAG_Y );  diag.push_back( DIAG_VY );
        diag.push_back( DIAG_Z );  diag.push_back( DIAG_VZ );
        diag.push_back( DIAG_CURR );
        TrajectoryDiagnosticData tdata;
        pdb.trajectories_at_plane( tdata, AXIS_Z, z_diag, diag );
        const size_t N = tdata.traj_count();

        double Itot = 0.0, sx=0, sx2=0, sy=0, sy2=0; size_t n=0;
        for( size_t i = 0; i < N; ++i ) {
            double vx = tdata(i,1).val(), vy = tdata(i,3).val(), vz = tdata(i,5).val();
            double I  = tdata(i,6).val();
            if( vz <= 0.0 ) continue;
            double xp = vx/vz, yp = vy/vz;               // divergence angles [rad]
            Itot += I; sx+=xp; sx2+=xp*xp; sy+=yp; sy2+=yp*yp; ++n;
        }
        double rms_x = n? std::sqrt(std::max(0.0, sx2/n - (sx/n)*(sx/n))) : 0.0;
        double rms_y = n? std::sqrt(std::max(0.0, sy2/n - (sy/n)*(sy/n))) : 0.0;
        double half_angle_deg = std::sqrt(rms_x*rms_x + rms_y*rms_y) * 180.0/M_PI;

        std::ofstream js( opath("diagnostics.json").c_str() );
        js << "{\n"
           << "  \"apertures_half\": " << aps.size() << ",\n"
           << "  \"z_diag_m\": " << z_diag << ",\n"
           << "  \"n_traj_at_plane\": " << n << ",\n"
           << "  \"I_at_plane_half_A\": " << Itot << ",\n"
           << "  \"I_full_grid_A\": " << 2.0*Itot << ",\n"
           << "  \"rms_div_x_rad\": " << rms_x << ",\n"
           << "  \"rms_div_y_rad\": " << rms_y << ",\n"
           << "  \"merged_half_angle_deg\": " << half_angle_deg << "\n"
           << "}\n";
        js.close();
        std::printf("diagnostics: I_full=%.4e A  merged half-angle=%.3f deg  (z=%.4f)\n",
                    2.0*Itot, half_angle_deg, z_diag);
    } catch( Error e ) {
        std::printf("WARNING: diagnostics block failed (reconcile API on cluster)\n");
        e.print_error_message( ibsimu.message(0) );
    }
}

int main( int argc, char **argv )
{
    try {
        simu( argc, argv );
    } catch( Error e ) {
        e.print_error_message( ibsimu.message(0) );
        return 1;
    }
    return 0;
}
