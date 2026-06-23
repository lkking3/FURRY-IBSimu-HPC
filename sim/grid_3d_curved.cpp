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
#include "particles.hpp"             // ParticleP3D (trajectory points)
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
    // Reference axis not parallel to n
    const double ax = ( std::fabs(n[0]) < 0.9 ) ? 1.0 : 0.0;
    const double ay = ( std::fabs(n[0]) < 0.9 ) ? 0.0 : 1.0;
    const double az = 0.0;
    // d1 = a x n  (then normalise)
    double d1x = ay*n[2] - az*n[1];
    double d1y = az*n[0] - ax*n[2];
    double d1z = ax*n[1] - ay*n[0];
    double l1 = std::sqrt( d1x*d1x + d1y*d1y + d1z*d1z );
    d1 = Vec3D( d1x/l1, d1y/l1, d1z/l1 );
    // d2 = n x d1  (already unit length: n, d1 orthonormal)  -> d1 x d2 = n
    double d2x = n[1]*d1[2] - n[2]*d1[1];
    double d2y = n[2]*d1[0] - n[0]*d1[2];
    double d2z = n[0]*d1[1] - n[1]*d1[0];
    double l2 = std::sqrt( d2x*d2x + d2y*d2y + d2z*d2z );
    d2 = Vec3D( d2x/l2, d2y/l2, d2z/l2 );
}

// Write a MeshScalarField as a legacy-VTK STRUCTURED_POINTS volume (BINARY,
// big-endian as the format requires). 'stride' downsamples to keep files small
// for ParaView (stride 2 -> 1/8 the nodes). x varies fastest.
static float to_be_float( double d )
{
    float v = (float)d; uint32_t x;
    std::memcpy( &x, &v, 4 ); x = __builtin_bswap32( x );
    std::memcpy( &v, &x, 4 ); return v;
}
static void write_vtk_scalar( const std::string &path, const char *name,
                              MeshScalarField &f, int s, bool mirror )
{
    if( s < 1 ) s = 1;
    const uint32_t nx = f.size(0), ny = f.size(1), nz = f.size(2);
    const double ox = f.origo(0), oy = f.origo(1), oz = f.origo(2), h = f.h();
    const uint32_t dx = (nx-1)/s + 1, dy = (ny-1)/s + 1, dz = (nz-1)/s + 1;
    // mirror across the x = ox plane (the x=0 symmetry boundary) -> full beam
    const uint32_t Dx = mirror ? (2*dx - 1) : dx;
    const double   Ox = mirror ? ( ox - (double)(dx-1)*s*h ) : ox;
    std::ofstream o( path.c_str(), std::ios::binary );
    o << "# vtk DataFile Version 3.0\n"
      << "grid_3d_curved " << name << "\nBINARY\n"
      << "DATASET STRUCTURED_POINTS\n"
      << "DIMENSIONS " << Dx << " " << dy << " " << dz << "\n"
      << "ORIGIN "  << Ox << " " << oy << " " << oz << "\n"
      << "SPACING " << h*s << " " << h*s << " " << h*s << "\n"
      << "POINT_DATA " << (size_t)Dx*dy*dz << "\n"
      << "SCALARS " << name << " float 1\nLOOKUP_TABLE default\n";
    for( uint32_t k = 0; k < nz; k += s )
        for( uint32_t j = 0; j < ny; j += s ) {
            if( mirror )
                for( int m = -(int)(dx-1); m <= (int)(dx-1); ++m ) {
                    uint32_t i = (uint32_t)( std::abs(m) * s );
                    float v = to_be_float( f(i,j,k) );
                    o.write( (char*)&v, 4 );
                }
            else
                for( uint32_t i = 0; i < nx; i += s ) {
                    float v = to_be_float( f(i,j,k) );
                    o.write( (char*)&v, 4 );
                }
        }
    o.close();
}

// Write saved particle trajectories as a legacy-VTK POLYDATA line set so the
// actual beamlets can be viewed in 3-D in ParaView (no display needed). Every
// 'stride'-th trajectory is written to keep the file manageable.
static void write_vtk_trajectories( const std::string &path,
                                    ParticleDataBase3D &pdb, int stride,
                                    bool mirror, double x0 )
{
    if( stride < 1 ) stride = 1;
    std::vector<uint32_t> sel;
    for( uint32_t i = 0; i < pdb.size(); i += stride )
        if( pdb.traj_size(i) >= 2 ) sel.push_back( i );
    size_t npts = 0;
    for( uint32_t i : sel ) npts += pdb.traj_size( i );
    const int blocks = mirror ? 2 : 1;  // block 1 = original, block 2 = x-mirror

    std::ofstream o( path.c_str() );
    o << "# vtk DataFile Version 3.0\nbeamlet trajectories\nASCII\n"
      << "DATASET POLYDATA\nPOINTS " << (npts*blocks) << " float\n";
    for( int b = 0; b < blocks; ++b )
        for( uint32_t i : sel ) {
            const uint32_t n = pdb.traj_size( i );
            for( uint32_t j = 0; j < n; ++j ) {
                const ParticleP3D &p =
                    static_cast<const ParticleP3D&>( pdb.trajectory_point( i, j ) );
                Vec3D L = p.location();
                double X = (b == 0) ? L[0] : (2.0*x0 - L[0]);
                o << X << " " << L[1] << " " << L[2] << "\n";
            }
        }
    const size_t nlines = sel.size()*blocks;
    o << "LINES " << nlines << " " << (npts*blocks + nlines) << "\n";
    size_t off = 0;
    for( int b = 0; b < blocks; ++b )
        for( uint32_t i : sel ) {
            const uint32_t n = pdb.traj_size( i );
            o << n;
            for( uint32_t j = 0; j < n; ++j ) o << " " << (off + j);
            o << "\n";
            off += n;
        }
    o.close();
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

    const std::string SCREEN_STL = envs("SCREEN_STL","geometry/screen_pos.stl");
    const std::string ACCEL_STL  = envs("ACCEL_STL", "geometry/accel_pos.stl");
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
    const double SC_ALPHA = envd("SC_ALPHA", 0.5);      // numerical under-relaxation
    // Physical space-charge compensation in the field-free drift (carried over
    // from the 2-D module's SC_FACTOR): downstream of SC_RAMP_START_Z the beam's
    // own charge is scaled by SC_FACTOR (1 = full SC, ->0 = fully neutralised),
    // ramped over SC_RAMP_LEN_Z (0 = step).
    const double SC_FACTOR       = envd("SC_FACTOR", 0.0005);     // match 2-D default
    const double SC_RAMP_START_Z = envd("SC_RAMP_START_Z", 0.050);// just past accel (m)
    const double SC_RAMP_LEN_Z   = envd("SC_RAMP_LEN_Z", 0.0);    // ramp length (m)

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
    const size_t ap_denom = std::max( (size_t)1, aps.size() );
    const int n_per = std::max( 1, (int)( (size_t)NPART / ap_denom ) );
    const double Te_J = TE * QE;
    const double mi   = M_AMU * AMU;
    const double cs   = std::sqrt( std::max(Te_J/mi, 0.0) );        // Bohm speed
    const double J    = JSC * 0.6 * Q_E * QE * NI * cs;             // [A/m^2]
    const double E0   = std::max( 1.0, std::fabs(VS - VA) );        // [eV]

    std::printf("apertures=%zu  n_per=%d  J=%.3e A/m^2  E0=%.1f eV  species: q=%g m=%g amu\n",
                aps.size(), n_per, J, E0, Q_E, M_AMU );

    // Downstream diagnostic plane + per-iteration convergence history.
    const double z_diag = envd("Z_DIAG", LZ - 5.0*H);
    std::vector<double> hist_I, hist_div; std::vector<size_t> hist_n;

    // Current-weighted transmitted current [A, half] and RMS divergence [rad] at
    // a z = const plane. Throws if no trajectory crosses the plane yet.
    auto plane_diag = [&]( double zpl, double &I_half, size_t &nplane, double &rms ){
        std::vector<trajectory_diagnostic_e> diag;
        diag.push_back( DIAG_XP ); diag.push_back( DIAG_YP ); diag.push_back( DIAG_CURR );
        TrajectoryDiagnosticData td;
        pdb.trajectories_at_plane( td, AXIS_Z, zpl, diag );
        const std::vector<double> &xp = td(0).data();
        const std::vector<double> &yp = td(1).data();
        const std::vector<double> &cu = td(2).data();
        nplane = cu.size(); I_half = 0.0; double s2 = 0.0;
        for( size_t i = 0; i < nplane; ++i ) { I_half += cu[i]; s2 += cu[i]*(xp[i]*xp[i]+yp[i]*yp[i]); }
        rms = (I_half > 0.0) ? std::sqrt( s2 / I_half ) : 0.0;
    };

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

        // per-iteration convergence readout at the downstream plane
        { double Ih=0.0, rmsd=0.0; size_t npl=0;
          try {
              plane_diag( z_diag, Ih, npl, rmsd );
              std::printf("[iter %d/%d]  I_full=%.4e A  n@plane=%zu  half-angle=%.3f deg\n",
                          it+1, ITERM, 2.0*Ih, npl, rmsd*180.0/M_PI);
              hist_I.push_back(2.0*Ih); hist_div.push_back(rmsd*180.0/M_PI); hist_n.push_back(npl);
          } catch( ... ) {
              std::printf("[iter %d/%d]  (diag skipped: no plane crossings yet)\n", it+1, ITERM);
          }
        }
        std::fflush(stdout);

        // physical space-charge compensation in the field-free drift (carried
        // over from 2-D): scale the beam charge by SC_FACTOR for z >= start.
        if( SC_FACTOR != 1.0 ) {
            const uint32_t nx = scharge.size(0), ny = scharge.size(1), nz = scharge.size(2);
            const double z0 = scharge.origo(2), hz = scharge.h();
            for( uint32_t k = 0; k < nz; ++k ) {
                const double z = z0 + hz*k;
                double f = 1.0;
                if( z >= SC_RAMP_START_Z ) {
                    if( SC_RAMP_LEN_Z > 0.0 ) {
                        double t = std::clamp( (z - SC_RAMP_START_Z)/SC_RAMP_LEN_Z, 0.0, 1.0 );
                        f = 1.0 + (SC_FACTOR - 1.0)*t;
                    } else f = SC_FACTOR;
                }
                if( f != 1.0 )
                    for( uint32_t j = 0; j < ny; ++j )
                        for( uint32_t i = 0; i < nx; ++i )
                            scharge(i,j,k) *= f;
            }
        }

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

    // ----------------------------- final diagnostics -----------------------
    // Headline numbers = last iteration's plane readout; convergence history is
    // written so you can see whether current/divergence had settled.
    {
        const double Itot  = hist_I.empty()   ? 0.0 : hist_I.back();   // full grid [A]
        const double hadeg = hist_div.empty() ? 0.0 : hist_div.back(); // merged half-angle [deg]
        const size_t Npl   = hist_n.empty()   ? 0   : hist_n.back();

        std::ofstream js( opath("diagnostics.json").c_str() );
        js << "{\n"
           << "  \"apertures_half\": " << aps.size() << ",\n"
           << "  \"z_diag_m\": " << z_diag << ",\n"
           << "  \"iterations\": " << hist_I.size() << ",\n"
           << "  \"n_traj_at_plane\": " << Npl << ",\n"
           << "  \"I_full_grid_A\": " << Itot << ",\n"
           << "  \"merged_half_angle_deg\": " << hadeg << ",\n"
           << "  \"history_I_full_A\": [";
        for( size_t i = 0; i < hist_I.size(); ++i )   js << (i?", ":"") << hist_I[i];
        js << "],\n  \"history_half_angle_deg\": [";
        for( size_t i = 0; i < hist_div.size(); ++i ) js << (i?", ":"") << hist_div[i];
        js << "],\n  \"history_n_at_plane\": [";
        for( size_t i = 0; i < hist_n.size(); ++i )   js << (i?", ":"") << hist_n[i];
        js << "]\n}\n";
        js.close();
        std::printf("FINAL: I_full=%.4e A  half-angle=%.3f deg  (z=%.4f, %zu iters)\n",
                    Itot, hadeg, z_diag, hist_I.size());
    }

    // ----------------------------- axial envelope scan ---------------------
    // Beam radius and divergence vs z, so the waist/focus is located and the
    // merged divergence is read at a meaningful plane (not an arbitrary one).
    if( envi("WRITE_ENVELOPE", 1) ) {
        try {
            const int    NZ = envi("ENV_NZ", 60);
            const double z0 = envd("ENV_Z0", 0.046);        // just downstream of accel
            const double z1 = envd("ENV_Z1", LZ - 2.0*H);
            std::ofstream es( opath("envelope.csv").c_str() );
            es << "z_m,I_full_A,r_rms_m,r95_m,r_max_m,half_angle_deg\n";
            for( int iz = 0; iz < NZ; ++iz ) {
                const double zz = z0 + (z1 - z0) * iz / std::max(1, NZ-1);
                std::vector<trajectory_diagnostic_e> dg;
                dg.push_back(DIAG_X);  dg.push_back(DIAG_Y);
                dg.push_back(DIAG_XP); dg.push_back(DIAG_YP); dg.push_back(DIAG_CURR);
                TrajectoryDiagnosticData td;
                try { pdb.trajectories_at_plane( td, AXIS_Z, zz, dg ); }
                catch( ... ) { continue; }
                const std::vector<double> &X  = td(0).data();
                const std::vector<double> &Y  = td(1).data();
                const std::vector<double> &XP = td(2).data();
                const std::vector<double> &YP = td(3).data();
                const std::vector<double> &CU = td(4).data();
                double I = 0.0, sr2 = 0.0, sdiv2 = 0.0;
                std::vector<double> rr;  rr.reserve( CU.size() );
                for( size_t i = 0; i < CU.size(); ++i ) {
                    const double r = std::sqrt( X[i]*X[i] + Y[i]*Y[i] );
                    const double w = CU[i];
                    I += w; sr2 += w*r*r; sdiv2 += w*(XP[i]*XP[i] + YP[i]*YP[i]);
                    rr.push_back(r);
                }
                double rrms = (I>0.0) ? std::sqrt(sr2/I) : 0.0;
                double hadeg= (I>0.0) ? std::sqrt(sdiv2/I)*180.0/M_PI : 0.0;
                double rmax = 0.0, r95 = 0.0;
                if( !rr.empty() ) { std::sort(rr.begin(), rr.end());
                    rmax = rr.back(); r95 = rr[(size_t)(0.95*(rr.size()-1))]; }
                es << zz << "," << 2.0*I << "," << rrms << "," << r95 << ","
                   << rmax << "," << hadeg << "\n";
            }
            es.close();
            std::printf("wrote envelope.csv (%d planes, z=%.3f..%.3f m)\n", NZ, z0, z1);
        } catch( Error e ) {
            std::printf("WARNING: envelope scan failed\n");
            e.print_error_message( ibsimu.message(0) );
        }
    }

    // ----------------------------- visualisation (headless PNG) ------------
    // GeomPlotter renders 2-D cut-planes of the 3-D model straight to PNG -- no
    // display/X11 needed. Data is already saved above, so a plotting failure is
    // non-fatal.
    if( envi("WRITE_PNG", 1) ) {
        try {
            const int W    = envi("PNG_W", 1800), Hpx = envi("PNG_H", 900);
            const int ylev = (int)( scharge.size(1) / 2 );          // y=0 slice (mirror plane)
            const int zlev = (int)std::lround( z_diag / H );        // z = z_diag slice

            // (a) longitudinal view: beamlets through the grids + drift
            GeomPlotter gp( geom );
            gp.set_size( W, Hpx );
            gp.set_view( VIEW_ZX, ylev );
            gp.set_epot( &epot );
            gp.set_efield( &efield );
            gp.set_particle_database( &pdb );
            gp.plot_png( opath("beam_zx.png").c_str() );

            // (b) close-up on the grid stack  (z in [0,0.06] m, x in [-0.05,0.05] m)
            GeomPlotter gz( geom );
            gz.set_size( 1200, 1000 );
            gz.set_view( VIEW_ZX, ylev );
            gz.set_epot( &epot );
            gz.set_efield( &efield );
            gz.set_particle_database( &pdb );
            gz.set_ranges( 0.0, -0.05, 0.06, 0.05 );
            gz.plot_png( opath("beam_zx_closeup.png").c_str() );

            // (c) transverse beam footprint at the diagnostic plane
            GeomPlotter gt( geom );
            gt.set_size( 1000, 1000 );
            gt.set_view( VIEW_XY, zlev );
            gt.set_epot( &epot );
            gt.set_particle_database( &pdb );
            gt.plot_png( opath("beam_xy_plane.png").c_str() );

            std::printf("wrote PNGs: beam_zx.png, beam_zx_closeup.png, beam_xy_plane.png\n");
        } catch( Error e ) {
            std::printf("WARNING: PNG plotting failed (state .dat files already saved)\n");
            e.print_error_message( ibsimu.message(0) );
        }
    }

    // ----------------------------- VTK volumes (ParaView 3-D) --------------
    // Potential + beam-density as STRUCTURED_POINTS volumes. Downsample with
    // VTK_STRIDE to keep files small (default 2 -> 1/8 nodes). Set WRITE_VTK=0
    // to skip (the density field allocates one extra mesh-sized array).
    if( envi("WRITE_VTK", 1) ) {
        try {
            const int  stride = envi("VTK_STRIDE", 2);
            const bool mirror = envi("VTK_MIRROR", 1) != 0;   // reflect half -> full beam
            write_vtk_scalar( opath("epot.vtk"), "potential", epot, stride, mirror );
            MeshScalarField tdens( geom );
            pdb.build_trajectory_density_field( tdens );
            write_vtk_scalar( opath("beam_density.vtk"), "beam_density", tdens, stride, mirror );
            // actual beamlet lines for a 3-D view in ParaView (mirror plane x=0)
            write_vtk_trajectories( opath("beam_trajectories.vtk"), pdb,
                                    envi("TRAJ_STRIDE", 4), mirror, 0.0 );
            std::printf("wrote VTK: epot.vtk, beam_density.vtk, beam_trajectories.vtk "
                        "(stride=%d, mirror=%d)\n", stride, (int)mirror);
        } catch( Error e ) {
            std::printf("WARNING: VTK export failed (state .dat files already saved)\n");
            e.print_error_message( ibsimu.message(0) );
        }
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
