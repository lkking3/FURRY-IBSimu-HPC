#include "geometry.hpp"
#include "func_solid.hpp"
#include "geomplotter.hpp"
#include "epot_bicgstabsolver.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "particledatabase.hpp"   // ParticleDataBase2D
#include "trajectorydiagnostics.hpp"
#include "config.h"               // for InitialPlasma
#include "meshscalarfield.hpp"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <string>
#include <ctime>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <iomanip>   // PATCH: for CSV precision
#include <vector>
#include <limits>
#include <chrono>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ---------- env helpers ----------
static double envd(const char* k, double vdef) {
    if (const char* s = std::getenv(k)) {
        char* e = nullptr; double v = std::strtod(s, &e);
        if (e && *e == '\0') return v;
    }
    return vdef;
}
static int envi(const char* k, int vdef) {
    if (const char* s = std::getenv(k)) return std::atoi(s);
    return vdef;
}
static inline double deg2rad(double deg){ return deg * (M_PI/180.0); }
static inline double rad2deg(double r){ return r * (180.0/M_PI); }

// v^(3/2) helper (avoids pow for speed + stability)
static inline double pow32(double v) { return (v > 0.0) ? (v * std::sqrt(v)) : 0.0; }

// Estimate number of beamlets (apertures) for a full grid from aperture radius (m).
// Table from user packing chart (radius_mm -> max apertures). Piecewise linear interpolation.
static double estimate_beamlet_count(double ap_rad_m) {
    struct Pt { double r; double n; };
    static const Pt T[] = {
        {1.0e-3, 2221.0}, {2.0e-3,  769.0}, {3.0e-3,  379.0}, {4.0e-3, 223.0},
        {5.0e-3,  151.0}, {6.0e-3,  109.0}, {7.0e-3,   73.0}, {8.0e-3,  55.0}
    };
    const int N = (int)(sizeof(T)/sizeof(T[0]));
    if (!(ap_rad_m > 0.0)) return T[0].n;
    if (ap_rad_m <= T[0].r) return T[0].n;
    if (ap_rad_m >= T[N-1].r) return T[N-1].n;
    for (int i = 0; i < N-1; ++i) {
        if (ap_rad_m >= T[i].r && ap_rad_m <= T[i+1].r) {
            const double r0=T[i].r, n0=T[i].n;
            const double r1=T[i+1].r, n1=T[i+1].n;
            const double tt = (r1>r0) ? ((ap_rad_m - r0)/(r1 - r0)) : 0.0;
            return n0 + tt*(n1 - n0);
        }
    }
    return T[N-1].n;
}


// env string helper
static std::string envs(const char* k, const char* defv = nullptr) {
    const char* s = std::getenv(k);
    if (s && *s) return std::string(s);
    return defv ? std::string(defv) : std::string();
}

// mkdir -p helpers
static int mkdir_one(const std::string &p) {
    if (p.empty()) return 0;
    if (::mkdir(p.c_str(), 0775) == 0) return 0;
    if (errno == EEXIST) return 0;
    return -1;
}
static int mkdir_p(const std::string &path) {
    if (path.empty()) return 0;
    std::string cur;
    if (path[0] == '/') { cur = "/"; }
    for (size_t i = (path[0]=='/'); i < path.size(); ++i) {
        char c = path[i];
        cur.push_back(c);
        if (c=='/') { if (mkdir_one(cur)) return -1; }
    }
    return mkdir_one(cur);
}

// timestamp + basename for output files
static std::string timestamp_now() {
    char buf[32];
    std::time_t t = std::time(nullptr);
    std::tm tmv{};
    #ifdef _WIN32
    localtime_s(&tmv, &t);
    #else
    tmv = *std::localtime(&t);
    #endif
    std::strftime(buf, sizeof(buf), "%Y%m%d_%H%M%S", &tmv);
    return std::string(buf);
}
static std::string basename_only(const std::string &p) {
    size_t pos = p.find_last_of('/');
    if (pos == std::string::npos) return p;
    return p.substr(pos+1);
}
// add suffix before file extension, preserving directory
static std::string add_suffix_before_ext(const std::string &path, const std::string &suffix) {
    size_t slash = path.find_last_of('/');
    std::string dir  = (slash == std::string::npos) ? "" : path.substr(0, slash+1);
    std::string file = (slash == std::string::npos) ? path : path.substr(slash+1);
    size_t dot = file.find_last_of('.');
    std::string name = (dot == std::string::npos) ? file : file.substr(0, dot);
    std::string ext  = (dot == std::string::npos) ? ""   : file.substr(dot);
    return dir + name + suffix + ext;
}



// ---------- globals for FuncSolid ----------
static double g_a        = 0.0;   // base aperture radius [m]
static double g_xs0      = 0.0;   // screen start x
static double g_xs1      = 0.0;   // screen end x
static double g_xa0      = 0.0;   // accel start x
static double g_xa1      = 0.0;   // accel end x
static double g_off      = 0.0;   // accelerator aperture y-offset [m]
static double g_scr_off  = 0.0;   // screen aperture y-offset [m]

// screen chamfers (upstream & downstream)
static double g_scr_du = 0.0, g_scr_mu = 0.0; // depth & slope upstream
static double g_scr_dd = 0.0, g_scr_md = 0.0; // depth & slope downstream
// accel chamfers (upstream & downstream)
static double g_acc_du = 0.0, g_acc_mu = 0.0;
static double g_acc_dd = 0.0, g_acc_md = 0.0;

// tube/endplate geometry globals
static double g_ybox = 0.0;
static double g_tube_x0 = 0.0, g_tube_x1 = 0.0;
static double g_tube_w  = 0.0;  // wall thickness
static double g_end_w   = 0.0;  // endplate thickness



// clamp two depths to fit within thickness (proportional if sum>t)
static inline void clamp_biface(double t, double &du, double &dd, bool &scaled) {
    du = std::clamp(du, 0.0, t);
    dd = std::clamp(dd, 0.0, t);
    if (du + dd > t && (du+dd)>0.0) {
        double f = t / (du + dd);
        du *= f; dd *= f;
        scaled = true;
    }
}

// half-height of open slot (aperture) as function of x inside each plate
static inline double screen_a_eff(double x) {
    double ae = g_a;
    // upstream face at x = g_xs0, chamfer zone: [xs0, xs0+du]
    if (g_scr_du > 0.0 && g_scr_mu > 0.0) {
        const double xu_end = g_xs0 + g_scr_du;
        if (x >= g_xs0 && x <= xu_end) ae += (xu_end - x) * g_scr_mu;
    }
    // downstream face at x = g_xs1, chamfer zone: [xs1-dd, xs1]
    if (g_scr_dd > 0.0 && g_scr_md > 0.0) {
        const double xd_sta = g_xs1 - g_scr_dd;
        if (x >= xd_sta && x <= g_xs1) ae += (x - xd_sta) * g_scr_md;
    }
    return ae;
}
static inline double accel_a_eff(double x) {
    double ae = g_a;
    // upstream face at x = g_xa0, chamfer zone: [xa0, xa0+du]
    if (g_acc_du > 0.0 && g_acc_mu > 0.0) {
        const double xu_end = g_xa0 + g_acc_du;
        if (x >= g_xa0 && x <= xu_end) ae += (xu_end - x) * g_acc_mu;
    }
    // downstream face at x = g_xa1, chamfer zone: [xa1-dd, xa1]
    if (g_acc_dd > 0.0 && g_acc_md > 0.0) {
        const double xd_sta = g_xa1 - g_acc_dd;
        if (x >= xd_sta && x <= g_xa1) ae += (x - xd_sta) * g_acc_md;
    }
    return ae;
}

// --------- solids: grids ----------
static bool screen_fn(double x, double y, double /*z*/) {
    if (x >= g_xs0 && x <= g_xs1) return (std::fabs(y - g_scr_off) >= screen_a_eff(x));
    return false;
}
static bool accel_fn(double x, double y, double /*z*/) {
    if (x >= g_xa0 && x <= g_xa1) return (std::fabs(y - g_off) >= accel_a_eff(x));
    return false;
}

// --------- solids: downstream tube walls (top & bottom) ----------
static bool tube_top_fn(double x, double y, double /*z*/) {
    return (x >= g_tube_x0 && x <= g_tube_x1 && y >= (g_ybox - g_tube_w));
}
static bool tube_bot_fn(double x, double y, double /*z*/) {
    return (x >= g_tube_x0 && x <= g_tube_x1 && y <= (-g_ybox + g_tube_w));
}

// --------- solid: downstream endplate (sample face) ----------
static bool endplate_fn(double x, double y, double /*z*/) {
    return (x >= (g_tube_x1 - g_end_w) && x <= g_tube_x1);
}

int main() try {

// ---- timing ----
using SteadyClock = std::chrono::steady_clock;
const auto t_start = SteadyClock::now();
auto t_after_mesh  = t_start;
auto t_after_solve = t_start;
auto t_after_diag  = t_start;
auto t_after_png   = t_start;
auto t_end         = t_start;

auto sec_since = [](const SteadyClock::time_point &a, const SteadyClock::time_point &b) -> double {
    return std::chrono::duration<double>(b - a).count();
};

    // ----- base geometry (env overrides) -----
    const double h        = envd("H",             1.0e-4);
    const double a        = envd("AP_RAD_M",      1.0e-3);
    const double t        = envd("GRID_T_M",      5.0e-4);
    const double gap      = envd("GAP_M",         3.0e-3);
    const double x_left   = envd("X_LEFT_M",      2.0e-3);
    const double x_right  = envd("X_RIGHT_M",     8.0e-3);
    const double x_right_phys = envd("X_RIGHT_PHYS_M", x_right);
    const bool truncated_drift = (x_right + 1.0e-12 < x_right_phys);
    const double off_y    = envd("ACCEL_OFF_Y_M", 0.0);
    const double scr_off_y = envd("SCREEN_OFF_Y_M", envd("PG_OFF_Y_M", 0.0));

    // ----- electrical potentials -----
    const double VS_V     = envd("VS_V", 0.0);   // screen potential [V]
    const double VA_V     = envd("VA_V",  -1.0e4);     // accel/tube potential [V]
    const double SAMPLE_V = envd("SAMPLE_V", VA_V); // sample face potential [V]

    // ----- screen chamfers -----
    double scr_du = envd("SCR_UP_DEPTH_M",   0.0);
    double scr_ua = envd("SCR_UP_ANGLE_DEG", 0.0);
    double scr_dd = envd("SCR_DN_DEPTH_M",   0.0);
    double scr_da = envd("SCR_DN_ANGLE_DEG", 0.0);

    // ----- accel chamfers -----
    double acc_du = envd("ACC_UP_DEPTH_M",   0.0);
    double acc_ua = envd("ACC_UP_ANGLE_DEG", 0.0);
    double acc_dd = envd("ACC_DN_DEPTH_M",   0.0);
    double acc_da = envd("ACC_DN_ANGLE_DEG", 0.0);

    // Toggle ion / plasma physics
    const int ENABLE_IONS        = envi("ENABLE_IONS", 0);   // 0 = old behavior
    const int ION_ITER_MAX       = envi("ION_ITER_MAX", 8);  // Vlasov/plasma iterations

    // Ion species (you can tune later)
    const double ION_Q_E         = envd("ION_Q_E",   1.0);        // charge in units of e
    const double ION_M_AMU       = envd("ION_M_AMU", 2.0);       // mass in atomic mass units (Ar+ default)
    const int    ION_NPART       = envi("ION_NPART", 2000);       // trajectories

    // Optional: transverse temperature of ions
    const double ION_TP_EV       = envd("ION_TP_EV", 0.0);        // parallel temp [eV]
    const double ION_TT_EV       = envd("ION_TT_EV", 0.2);        // transverse temp [eV]

    // Plate positions
    const double xs0 = 0.0;
    const double xs1 = xs0 + t;
    const double xa0 = xs1 + gap;
    const double xa1 = xa0 + t;

    // ----- Perveance bookkeeping (Child-Langmuir reference + later geometric perveance) -----
    // Accelerating voltage magnitude used for perveance (V = |VS - VA|).
    const double V_accel_V = std::fabs(VS_V - VA_V);
    const double V32 = pow32(V_accel_V);

    // Effective gap used for CL estimate. Default = physical grid gap (screen dn face -> accel up face).
    // Override via PERV_D_EFF_M if you want a thickness/chamfer-corrected effective distance.
    const double d_eff_m = envd("PERV_D_EFF_M", gap);

    // Child-Langmuir (planar) current density: J_CL = (4/9) eps0 sqrt(2 q / m) * V^(3/2) / d^2
    // Perveance per aperture: P_CL = I_CL / V^(3/2) = (4/9) eps0 sqrt(2 q / m) * A / d^2
    const double eps0 = 8.8541878128e-12;      // [F/m]
    const double qe   = 1.602176634e-19;       // [C]
    const double amu  = 1.66053906660e-27;     // [kg]
    const double q_C  = std::fabs(ION_Q_E) * qe;
    const double m_kg = std::fabs(ION_M_AMU) * amu;
    const double A_ap_m2 = M_PI * a * a;

    double P_CL_A_per_V32 = 0.0;
    double I_CL_A = 0.0;
    if (d_eff_m > 0.0 && m_kg > 0.0 && q_C > 0.0) {
        P_CL_A_per_V32 = (4.0/9.0) * eps0 * std::sqrt(2.0*q_C/m_kg) * A_ap_m2 / (d_eff_m*d_eff_m);
        I_CL_A = P_CL_A_per_V32 * V32;
    }

    const double N_beamlets_est = estimate_beamlet_count(a);
    const double P_CL_sys_A_per_V32 = P_CL_A_per_V32 * N_beamlets_est;
    const double I_CL_sys_A = I_CL_A * N_beamlets_est;

    // Record key plasma/beam tuning knobs in meta.json even when ions are disabled (for provenance).
    const double PLASMA_NI_M3_META = envd("PLASMA_NI_M3", 5.0e18);
    const double PLASMA_TE_EV_META = envd("PLASMA_TE_EV", 3.7);
    const double PLASMA_UP_V_META  = envd("PLASMA_UP_V", VS_V);
    const double ION_J_SCALE_META     = envd("ION_J_SCALE", 1.0);
    const double SC_FACTOR_META       = envd("SC_FACTOR", 1.0);
    const double SC_RAMP_START_M_META = envd("SC_RAMP_START_M", xa1 + 0.5*h);
    const double SC_RAMP_LEN_M_META   = envd("SC_RAMP_LEN_M", 0.0);



    // Clamp chamfers to thickness
    bool scr_scaled=false, acc_scaled=false;
    clamp_biface(t, scr_du, scr_dd, scr_scaled);
    clamp_biface(t, acc_du, acc_dd, acc_scaled);

    // Convert angles to slopes
    const double scr_mu = std::tan(deg2rad(std::max(scr_ua, 0.0)));
    const double scr_md = std::tan(deg2rad(std::max(scr_da, 0.0)));
    const double acc_mu = std::tan(deg2rad(std::max(acc_ua, 0.0)));
    const double acc_md = std::tan(deg2rad(std::max(acc_da, 0.0)));

    // Push globals
    g_a = a;
    g_xs0 = xs0; g_xs1 = xs1;
    g_xa0 = xa0; g_xa1 = xa1;
    g_off = off_y;
    g_scr_off = scr_off_y;

    g_scr_du = scr_du; g_scr_mu = scr_mu;
    g_scr_dd = scr_dd; g_scr_md = scr_md;
    g_acc_du = acc_du; g_acc_mu = acc_mu;
    g_acc_dd = acc_dd; g_acc_md = acc_md;

    // Tube/endplate sizing
    const double max_delta = std::max({ scr_du*scr_mu, scr_dd*scr_md, acc_du*acc_mu, acc_dd*acc_md });
    const double off_abs   = std::max(std::fabs(off_y), std::fabs(scr_off_y));
    const double ybox_def  = std::max(3.0*(a + max_delta) + off_abs, 5.0e-3);
    g_ybox = envd("YBOX_M", ybox_def);

    // Domain extents
    const double xmin = -x_left;
    const double xmax =  xa1 + x_right;
    const double xmax_phys = xa1 + x_right_phys;

    // Downstream tube region spans from accel end to domain end
    g_tube_x0 = envd("TUBE_X_START_M", xa1);
    g_tube_x1 = envd("TUBE_X_END_M",   xmax);
    g_tube_w  = std::max(envd("TUBE_WALL_T_M", 2.0*h), h);    // keep >= h
    g_end_w   = std::max(envd("SAMPLE_PLATE_T_M", 2.0*h), h); // endplate thickness
    const bool endplate_enabled = !truncated_drift;
    
    // results directory naming (sweep-friendly)
    const std::string results_base = envs("RESULTS_DIR", "results");
    std::string stamp = envs("RUN_STAMP"); if (stamp.empty()) stamp = timestamp_now();
    const std::string prefix = envs("RUN_PREFIX", "run");
    const std::string tag    = envs("RUN_TAG", "");
    const std::string jobid  = envs("SLURM_JOB_ID", "");
    std::string runname = prefix + "_" + stamp;
    if (!tag.empty())   runname += "_" + tag;
    if (!jobid.empty()) runname += "_j" + jobid;
    const std::string outdir = results_base + "/" + runname;

    if (mkdir_p(results_base) || mkdir_p(outdir)) {
        std::fprintf(stderr, "ERROR: could not create results dir: %s\n", outdir.c_str());
        return 2;
    }

    // keep original PNG_NAME but write into results folder
    const char* pngname_c = std::getenv("PNG_NAME") ? std::getenv("PNG_NAME") : "two_grid_geom.png";
    const std::string png = outdir + "/" + basename_only(pngname_c);
    


    // Mesh size
    Int3D size( int(std::floor((xmax - xmin)/h)+1),
                int(std::floor((2.0*g_ybox)/h)+1),
                1 );

    Geometry geom( MODE_2D, size, Vec3D(xmin, -g_ybox, 0.0), h );

    // Define solids
    geom.set_solid(7,  new FuncSolid(&screen_fn));
    geom.set_solid(8,  new FuncSolid(&accel_fn));
    geom.set_solid(9,  new FuncSolid(&tube_top_fn));
    geom.set_solid(10, new FuncSolid(&tube_bot_fn));
    if (endplate_enabled) {
        geom.set_solid(11, new FuncSolid(&endplate_fn));
    }

    // --- Boundary conditions ---
    // Global box:
    geom.set_boundary(1, Bound(BOUND_DIRICHLET, VS_V)); // xmin (plasma chamber potential)
    if (endplate_enabled) {
        // Full-length run: pin far boundary to tube potential.
        geom.set_boundary(2, Bound(BOUND_DIRICHLET, VA_V)); // xmax
    } else {
        // Truncated drift: open boundary (no normal E-field).
        geom.set_boundary(2, Bound(BOUND_NEUMANN, 0.0)); // xmax
    }
    geom.set_boundary(3, Bound(BOUND_DIRICHLET,   VA_V));  // ymin (open upstream; tube walls provide BC downstream)
    geom.set_boundary(4, Bound(BOUND_DIRICHLET,   VA_V));  // ymax

    // Conductors:
    geom.set_boundary(7,  Bound(BOUND_DIRICHLET, VS_V));     // screen
    geom.set_boundary(8,  Bound(BOUND_DIRICHLET, VA_V));     // accel
    geom.set_boundary(9,  Bound(BOUND_DIRICHLET, VA_V));     // tube top
    geom.set_boundary(10, Bound(BOUND_DIRICHLET, VA_V));     // tube bottom
    if (endplate_enabled) {
        geom.set_boundary(11, Bound(BOUND_DIRICHLET, SAMPLE_V)); // endplate / sample face
    }

    geom.build_mesh();
    t_after_mesh = SteadyClock::now();


    // meta.json for sweep bookkeeping
    {
        std::ofstream m(outdir + "/meta.json");
        m << "{\n";
        // Provenance
        m << "  \"RUN_STAMP\":\"" << stamp << "\",\n";
        m << "  \"RUN_TAG\":\""   << envs("RUN_TAG","") << "\",\n";
        m << "  \"runname\":\""   << runname << "\",\n";

        // Geometry & domain
        m << "  \"H\":" << h << ",\n";
        m << "  \"AP_RAD_M\":" << a << ",\n";
        m << "  \"GRID_T_M\":" << t << ",\n";
        m << "  \"GAP_M\":" << gap << ",\n";
        m << "  \"ACCEL_OFF_Y_M\":" << off_y << ",\n";
        m << "  \"SCREEN_OFF_Y_M\":" << scr_off_y << ",\n";
        m << "  \"X_LEFT_M\":" << x_left << ",\n";
        m << "  \"X_RIGHT_M\":" << x_right << ",\n";
        m << "  \"X_RIGHT_PHYS_M\":" << x_right_phys << ",\n";
        m << "  \"TRUNCATED_DRIFT\":" << (truncated_drift ? "true" : "false") << ",\n";
        m << "  \"XMAX_PHYS_M\":" << xmax_phys << ",\n";
        m << "  \"YBOX_M\":" << g_ybox << ",\n";

        // Potentials (top-level + nested, for convenience)
        m << "  \"VS_V\":" << VS_V << ",\n";
        m << "  \"VA_V\":" << VA_V << ",\n";
        m << "  \"SAMPLE_V\":" << SAMPLE_V << ",\n";
        m << "  \"potentials\": {\"VS_V\":" << VS_V
        << ", \"VA_V\":" << VA_V << ", \"SAMPLE_V\":" << SAMPLE_V << "},\n";
        // Plasma/beam knobs (for provenance; used by Stage-2)
        m << "  \"physics\": {"
          << "\"ENABLE_IONS\":" << ENABLE_IONS
          << ", \"ION_ITER_MAX\":" << ION_ITER_MAX
          << ", \"ION_Q_E\":" << ION_Q_E
          << ", \"ION_M_AMU\":" << ION_M_AMU
          << ", \"ION_NPART\":" << ION_NPART
          << ", \"ION_TP_EV\":" << ION_TP_EV
          << ", \"ION_TT_EV\":" << ION_TT_EV
          << ", \"PLASMA_NI_M3\":" << PLASMA_NI_M3_META
          << ", \"PLASMA_TE_EV\":" << PLASMA_TE_EV_META
          << ", \"PLASMA_UP_V\":" << PLASMA_UP_V_META
          << ", \"ION_J_SCALE\":" << ION_J_SCALE_META
          << ", \"SC_FACTOR\":" << SC_FACTOR_META
          << ", \"SC_RAMP_START_M\":" << SC_RAMP_START_M_META
          << ", \"SC_RAMP_LEN_M\":" << SC_RAMP_LEN_M_META
          << "},\n";


        // Chamfers
        m << "  \"screen_chamfer\": {"
        << "\"up_depth\":" << scr_du << ", \"up_angle_deg\":" << scr_ua
        << ", \"dn_depth\":" << scr_dd << ", \"dn_angle_deg\":" << scr_da << "},\n";
        m << "  \"accel_chamfer\": {"
        << "\"up_depth\":" << acc_du << ", \"up_angle_deg\":" << acc_ua
        << ", \"dn_depth\":" << acc_dd << ", \"dn_angle_deg\":" << acc_da << "},\n";

        // Tube / endplate
        m << "  \"tube\": {\"x_start\":" << g_tube_x0 << ", \"x_end\":" << g_tube_x1
        << ", \"wall_t\":" << g_tube_w << ", \"endplate_t\":" << g_end_w << "},\n";

        // Perveance reference numbers (Child-Langmuir) for this geometry
        m << "  \"perveance_ref\": {"
          << "\"V_accel_V\":" << V_accel_V
          << ", \"V32\":" << V32
          << ", \"d_eff_m\":" << d_eff_m
          << ", \"A_ap_m2\":" << A_ap_m2
          << ", \"q_C\":" << q_C
          << ", \"m_kg\":" << m_kg
          << ", \"P_CL_A_per_V32\":" << P_CL_A_per_V32
          << ", \"I_CL_A\":" << I_CL_A
          << ", \"N_beamlets_est\":" << N_beamlets_est
          << ", \"P_CL_sys_A_per_V32\":" << P_CL_sys_A_per_V32
          << ", \"I_CL_sys_A\":" << I_CL_sys_A
          << "},\n";


        // Domain box
        m << "  \"domain\": {\"xmin\":" << xmin << ", \"xmax\":" << xmax
        << ", \"ymin\":" << -g_ybox << ", \"ymax\":" << g_ybox << "}\n";

        m << "}\n";
    }
    


    // Summary
    std::printf("# two-grid 2D (offset + chamfers + BCs + tube/endplate)\n");
    std::printf("VS_V=%g  VA_V=%g  SAMPLE_V=%g\n", VS_V, VA_V, SAMPLE_V);
    std::printf("screen: x=[%g,%g]  up(depth=%g,ang=%g)  dn(depth=%g,ang=%g)\n",
                xs0,xs1,scr_du,scr_ua,scr_dd,scr_da);
    std::printf("accel : x=[%g,%g]  up(depth=%g,ang=%g)  dn(depth=%g,ang=%g)\n",
                xa0,xa1,acc_du,acc_ua,acc_dd,acc_da);
    if (scr_scaled) std::printf("NOTE: screen chamfers scaled to fit thickness\n");
    if (acc_scaled) std::printf("NOTE: accel  chamfers scaled to fit thickness\n");
    std::printf("tube walls: x=[%g,%g]  wall_t=%g  endplate_t=%g\n",
                g_tube_x0,g_tube_x1,g_tube_w,g_end_w);
    std::printf("domain: x=[%g,%g], y=[%g,%g]\n", xmin,xmax,-g_ybox,g_ybox);
    std::printf("# results dir: %s\n", outdir.c_str());

    // Optional PNG of geometry (and epot contours if solved)
    const int writepng = envi("WRITE_PNG", 1);
    const char* pngname  = std::getenv("PNG_NAME") ? std::getenv("PNG_NAME") : "two_grid_geom.png";

    // Optional Laplace / plasma solve (optionally with ion beam)
    const int RUN_SOLVE   = envi("RUN_SOLVE", 0);

    if (RUN_SOLVE) {
        // Common fields
        EpotField       epot(geom);
        MeshScalarField scharge(geom);
        MeshVectorField bfield;         // required by iterate_trajectories()
        EpotEfield      efield(epot);
        EpotBiCGSTABSolver solver(geom);

        // Particle DB (only really used when ENABLE_IONS != 0)
        ParticleDataBase2D pdb(geom);
        bool use_pdb = false;

        if (!ENABLE_IONS) {
            // ----------------------------
            // Original behavior: Laplace
            // ----------------------------
            solver.solve(epot, scharge);
            efield.recalculate();

        } else {
            // ----------------------------------------
            // Plasma + ion beam with nonlinear plasma
            // ----------------------------------------

            // === User-requested plasma parameters ===
            const double PLASMA_NI_M3 = envd("PLASMA_NI_M3", 5.0e18);  // [m^-3]
            const double PLASMA_TE_EV = envd("PLASMA_TE_EV", 3.7);     // [eV]
            const double PLASMA_UP_V  = envd("PLASMA_UP_V", VS_V);     // [V] plasma potential (usually screen/source)

            // Ion species + beam parameters (override via env if needed)
            const double ION_Q_E      = envd("ION_Q_E",   1.0);        // charge in units of e
            const double ION_M_AMU    = envd("ION_M_AMU", 40.0);       // mass in amu (Ar+ default)
            const int    ION_NPART    = envi("ION_NPART", 2000);       // number of trajectories
            const int    ION_ITER_MAX = envi("ION_ITER_MAX", 8);       // Vlasov iterations

            const double ION_TP_EV    = envd("ION_TP_EV", 0.0);        // parallel temperature [eV]
            const double ION_TT_EV    = envd("ION_TT_EV", 0.2);        // transverse temperature [eV]


            // Perveance & neutralization
            const double ION_J_SCALE  = envd("ION_J_SCALE", 1.0);      // scale injected current
            const double SC_FACTOR    = envd("SC_FACTOR",   1.0);      // scale space charge (1 = full SC)
            const double SC_RAMP_START_M = SC_RAMP_START_M_META;       // start of SCC ramp (m)
            const double SC_RAMP_LEN_M   = SC_RAMP_LEN_M_META;         // ramp length (m), 0 = step

            use_pdb = true;

            // --- Nonlinear plasma model setup ---

            // Plasma volume: x < xs0 (upstream side of the screen)
            InitialPlasma initp(AXIS_X, xs0);
            solver.set_initial_plasma(PLASMA_UP_V, &initp);

            // Background electron density rho_e = -e n_e, used by positive exponential plasma model
            const double qe  = 1.602176634e-19;       // [C]
            const double kB  = 1.380649e-23;          // [J/K]
            const double amu = 1.66053906660e-27;     // [kg]

            const double ne   = PLASMA_NI_M3;         // assume quasi-neutral n_e ≈ n_i
            const double rhoe = -qe * ne;             // [C/m^3]

            // Enable nonlinear plasma
            solver.set_pexp_plasma(rhoe, PLASMA_TE_EV, PLASMA_UP_V);

            // --- Derive ion current density J from n_i, T_e (Bohm-like flux) ---
            const double Te_J  = PLASMA_TE_EV * qe;        // [J]
            const double mi_kg = ION_M_AMU * amu;
            double cs = 0.0;
            if (mi_kg > 0.0)
                cs = std::sqrt(std::max(Te_J / mi_kg, 0.0));   // ion sound speed

            // J ≈ 0.6 * q * n_i * c_s  [A/m^2], with user scale
            const double J_Apm2 = ION_J_SCALE * 0.6 * ION_Q_E * qe * ne * cs;


            // Full-domain model: no mirror boundaries
            bool pmirror[6] = { false, false, false, false, false, false };
            pdb.set_mirror(pmirror);

            // Beam start line: just upstream of the screen aperture, across its open height
            const double x_emit = std::max(xmin + 2.0*h, xs0 - 2.0*h);
            const double y_lo   = scr_off_y - a;
            const double y_hi   = scr_off_y + a;

            // Injection energy: potential drop between plasma and accel side
            const double E0_ev  = std::max(1.0, std::fabs(PLASMA_UP_V - VA_V));

            for (int it = 0; it < ION_ITER_MAX; ++it) {
                solver.solve(epot, scharge);
                efield.recalculate();

                pdb.clear();
                pdb.add_2d_beam_with_energy(
                    (uint32_t)ION_NPART,
                    J_Apm2,        // current density [A/m^2]
                    ION_Q_E,
                    ION_M_AMU,
                    E0_ev,         // beam energy [eV]
                    ION_TP_EV,     // parallel temperature [eV]
                    ION_TT_EV,     // transverse temperature [eV]
                    x_emit, y_lo,  // start of emission line
                    x_emit, y_hi   // end of emission line
                );
                pdb.iterate_trajectories(scharge, efield, bfield);

                // Emulate partial neutralization in the drift
                if (SC_FACTOR != 1.0) {
                    const uint32_t nx = scharge.size(0);
                    const uint32_t ny = scharge.size(1);
                    const double x0_sc = SC_RAMP_START_M;
                    const double L_sc = SC_RAMP_LEN_M;
                    const double inv_L = (L_sc > 0.0) ? (1.0 / L_sc) : 0.0;
                    const double x_origin = scharge.origo(0);
                    const double hx = scharge.h();

                    for (uint32_t i = 0; i < nx; ++i) {
                        const double x = x_origin + hx * i;
                        double f = 1.0;

                        if (x >= x0_sc) {
                            if (L_sc > 0.0) {
                                double t = (x - x0_sc) * inv_L;
                                t = std::clamp(t, 0.0, 1.0);
                                f = 1.0 + (SC_FACTOR - 1.0) * t;
                            } else {
                                f = SC_FACTOR;
                            }
                        }

                        if (f != 1.0) {
                            for (uint32_t j = 0; j < ny; ++j) {
                                scharge(i, j) *= f;
                            }
                        }
                    }
                }
            }
            // ====== Beam diagnostics: currents + collimation ======
            // Only meaningful when ions are enabled and pdb is populated.
            const double Leff = 0.785 * (2.0 * a);
            // Helper: total current and basic width at a given x-plane

            auto current_and_width_at_plane = [&](double xpos,
                                                double &Itot,
                                                double &y_mean,
                                                double &y_rms0,
                                                double &y_rms_c,
                                                double &y_absmax) {
                TrajectoryDiagnosticData tdata;
                std::vector<trajectory_diagnostic_e> diag;
                diag.push_back(DIAG_Y);
                diag.push_back(DIAG_CURR);

                // Collect diagnostic data at plane x = xpos
                pdb.trajectories_at_plane(tdata, AXIS_X, xpos, diag);

                const std::vector<double> &y    = tdata(0).data();
                const std::vector<double> &curr = tdata(1).data();

                Itot     = 0.0;
                y_mean   = 0.0;
                y_rms0   = 0.0;
                y_rms_c  = 0.0;
                y_absmax = 0.0;

                const size_t N = curr.size();
                double sum_wy = 0.0;
                double sum_wy2 = 0.0;
                for (size_t i = 0; i < N; ++i) {
                    const double w  = curr[i];      // current weight (A/m in 2D)
                    const double yy = y[i];
                    const double ay = std::fabs(yy);

                    Itot   += w;
                    sum_wy += w * yy;
                    sum_wy2 += w * yy * yy;
                    if (ay > y_absmax) y_absmax = ay;
                }

                if (Itot > 0.0) {
                    y_mean = sum_wy / Itot;
                    y_rms0 = std::sqrt(sum_wy2 / Itot);
                    const double var = (sum_wy2 / Itot) - (y_mean * y_mean);
                    y_rms_c = (var > 0.0) ? std::sqrt(var) : 0.0;
                } else {
                    y_mean = 0.0;
                    y_rms0 = 0.0;
                    y_rms_c = 0.0;
                }
            };

            // ---- 2.1 Currents at PG entrance, AG exit, sample wall ----

            // Planes:
            //  - PG entrance: upstream face of plasma/screen plate
            //  - AG exit:     downstream face of accel plate
            //  - Sample wall: just upstream of endplate
            // Probe planes (offset by +0.5*h so we sample *in vacuum*, not on an electrode surface)
            const double x_pg_plane = xs0 + 0.5*h;
            const double x_ag_plane = xa1 + 0.5*h;
            const double x_sm_plane = endplate_enabled ? (g_tube_x1 - g_end_w - 0.5*h)
                                                     : std::numeric_limits<double>::quiet_NaN();

            double I_pg_in_Apm  = 0.0, ymean_pg = 0.0, yr_pg  = 0.0, yr_pg_c  = 0.0, ymax_pg  = 0.0;
            double I_ag_out_Apm = 0.0, ymean_ag = 0.0, yr_ag  = 0.0, yr_ag_c  = 0.0, ymax_ag  = 0.0;
            double I_sm_Apm     = 0.0, ymean_sm = 0.0, yr_sm  = 0.0, yr_sm_c  = 0.0, ymax_sm  = 0.0;

            current_and_width_at_plane(x_pg_plane, I_pg_in_Apm,  ymean_pg,  yr_pg,  yr_pg_c,  ymax_pg);
            current_and_width_at_plane(x_ag_plane, I_ag_out_Apm, ymean_ag, yr_ag,  yr_ag_c,  ymax_ag);
            if (endplate_enabled) {

                current_and_width_at_plane(x_sm_plane, I_sm_Apm, ymean_sm, yr_sm, yr_sm_c, ymax_sm);

            }
            // 3D estimates in amps using Leff
            const double I_pg_in_A  = I_pg_in_Apm  * Leff;
            const double I_ag_out_A = I_ag_out_Apm * Leff;
            const double I_sm_A     = I_sm_Apm     * Leff;

            std::printf("[beam] currents: "
                        "PG entrance = %.3e A/m (%.3e A), "
                        "AG exit = %.3e A/m (%.3e A), "
                        "sample = %.3e A/m (%.3e A)\n",
                        I_pg_in_Apm,  I_pg_in_A,
                        I_ag_out_Apm, I_ag_out_A,
                        I_sm_Apm,     I_sm_A);

            // ---- 2.1a Sample radial/diameter profile ----
            {
                const double sample_diam_m = envd("SAMPLE_DIAM_M", 6.35e-2);
                const double sample_rad_m = 0.5 * sample_diam_m;
                const double sample_bin_m = envd("SAMPLE_BIN_M", 5.0e-4);
                const int allow_trunc_profile = envi("SAMPLE_PROFILE_ALLOW_TRUNC", 1);
                const double x_profile_env = envd("SAMPLE_PROFILE_X_M", std::numeric_limits<double>::quiet_NaN());
                const double x_profile_plane = std::isfinite(x_profile_env)
                                                   ? x_profile_env
                                                   : (endplate_enabled ? x_sm_plane : (xmax - 0.5*h));
                const bool profile_ok = endplate_enabled || (allow_trunc_profile != 0);

                if (profile_ok && sample_rad_m > 0.0 && sample_bin_m > 0.0 &&
                    x_profile_plane > (xmin + 0.5*h) && x_profile_plane < (xmax - 0.5*h)) {
                    TrajectoryDiagnosticData tdata;
                    std::vector<trajectory_diagnostic_e> diag;
                    diag.push_back(DIAG_Y);
                    diag.push_back(DIAG_CURR);
                    pdb.trajectories_at_plane(tdata, AXIS_X, x_profile_plane, diag);

                    const std::vector<double> &y    = tdata(0).data();
                    const std::vector<double> &curr = tdata(1).data();
                    const size_t N = curr.size();

                    const int nbins_r = (int)std::ceil(sample_rad_m / sample_bin_m);
                    const int nbins_y = (int)std::ceil((2.0 * sample_rad_m) / sample_bin_m);
                    if (nbins_r > 0 && nbins_y > 0) {
                        std::vector<double> bin_curr_r(nbins_r, 0.0);
                        std::vector<size_t> bin_cnt_r(nbins_r, 0);
                        std::vector<double> bin_curr_y(nbins_y, 0.0);
                        std::vector<size_t> bin_cnt_y(nbins_y, 0);
                        double sum_curr = 0.0;

                        for (size_t i = 0; i < N; ++i) {
                            const double yy = y[i];
                            const double r = std::fabs(yy);
                            if (r > sample_rad_m) {
                                continue;
                            }
                            const double w = curr[i];
                            sum_curr += w;

                            int idx_r = (int)std::floor(r / sample_bin_m);
                            if (idx_r < 0) idx_r = 0;
                            if (idx_r >= nbins_r) idx_r = nbins_r - 1;
                            bin_curr_r[idx_r] += w;
                            bin_cnt_r[idx_r] += 1;

                            const double y0 = -sample_rad_m;
                            int idx_y = (int)std::floor((yy - y0) / sample_bin_m);
                            if (idx_y < 0) idx_y = 0;
                            if (idx_y >= nbins_y) idx_y = nbins_y - 1;
                            bin_curr_y[idx_y] += w;
                            bin_cnt_y[idx_y] += 1;
                        }

                        std::ofstream pj(outdir + "/sample_radial_profile.json");
                        pj << "{\n";
                        pj << "  \"x_plane_m\": " << std::setprecision(10) << x_profile_plane << ",\n";
                        pj << "  \"sample_diam_m\": " << std::setprecision(10) << sample_diam_m << ",\n";
                        pj << "  \"sample_rad_m\": " << std::setprecision(10) << sample_rad_m << ",\n";
                        pj << "  \"bin_width_m\": " << std::setprecision(10) << sample_bin_m << ",\n";
                        pj << "  \"total_I_Apm\": " << std::setprecision(10) << sum_curr << ",\n";
                        pj << "  \"total_I_A\": " << std::setprecision(10) << (sum_curr * Leff) << ",\n";
                        pj << "  \"num_samples\": " << N << ",\n";
                        pj << "  \"bins\": [\n";
                        for (int b = 0; b < nbins_r; ++b) {
                            const double r0 = b * sample_bin_m;
                            const double r1 = (b + 1) * sample_bin_m;
                            pj << "    {\"r_lo_m\": " << std::setprecision(10) << r0
                               << ", \"r_hi_m\": " << std::setprecision(10) << r1
                               << ", \"I_Apm\": " << std::setprecision(10) << bin_curr_r[b]
                               << ", \"I_A\": " << std::setprecision(10) << (bin_curr_r[b] * Leff)
                               << ", \"fraction\": ";
                            if (sum_curr > 0.0) pj << std::setprecision(10) << (bin_curr_r[b] / sum_curr); else pj << "0";
                            pj << ", \"count\": " << bin_cnt_r[b] << "}";
                            if (b != nbins_r - 1) pj << ",";
                            pj << "\n";
                        }
                        pj << "  ]\n";
                        pj << "}\n";

                        std::ofstream dj(outdir + "/sample_diameter_profile.json");
                        dj << "{\n";
                        dj << "  \"x_plane_m\": " << std::setprecision(10) << x_profile_plane << ",\n";
                        dj << "  \"sample_diam_m\": " << std::setprecision(10) << sample_diam_m << ",\n";
                        dj << "  \"sample_rad_m\": " << std::setprecision(10) << sample_rad_m << ",\n";
                        dj << "  \"bin_width_m\": " << std::setprecision(10) << sample_bin_m << ",\n";
                        dj << "  \"total_I_Apm\": " << std::setprecision(10) << sum_curr << ",\n";
                        dj << "  \"total_I_A\": " << std::setprecision(10) << (sum_curr * Leff) << ",\n";
                        dj << "  \"num_samples\": " << N << ",\n";
                        dj << "  \"bins\": [\n";
                        for (int b = 0; b < nbins_y; ++b) {
                            const double y0 = -sample_rad_m + b * sample_bin_m;
                            const double y1 = y0 + sample_bin_m;
                            dj << "    {\"y_lo_m\": " << std::setprecision(10) << y0
                               << ", \"y_hi_m\": " << std::setprecision(10) << y1
                               << ", \"I_Apm\": " << std::setprecision(10) << bin_curr_y[b]
                               << ", \"I_A\": " << std::setprecision(10) << (bin_curr_y[b] * Leff)
                               << ", \"fraction\": ";
                            if (sum_curr > 0.0) dj << std::setprecision(10) << (bin_curr_y[b] / sum_curr); else dj << "0";
                            dj << ", \"count\": " << bin_cnt_y[b] << "}";
                            if (b != nbins_y - 1) dj << ",";
                            dj << "\n";
                        }
                        dj << "  ]\n";
                        dj << "}\n";
                    }
                } else if (profile_ok &&
                           (x_profile_plane <= (xmin + 0.5*h) || x_profile_plane >= (xmax - 0.5*h))) {
                    std::fprintf(stderr,
                                 "[beam] WARNING: sample profile plane out of bounds (x=%.6e m); skipping profile.\n",
                                 x_profile_plane);
                }
            }

            // Deflection (centroid steering) between AG exit and right boundary of computed space.
            const double x_right_plane = endplate_enabled
                                             ? (g_tube_x1 - g_end_w - 0.5*h)
                                             : (xmax - 0.5*h);
            double I_right_Apm = 0.0, ymean_right = 0.0, yr_right = 0.0, yr_right_c = 0.0, ymax_right = 0.0;
            bool have_right_plane = false;
            if (x_right_plane > x_ag_plane + 2.0*h) {
                current_and_width_at_plane(x_right_plane, I_right_Apm, ymean_right, yr_right, yr_right_c, ymax_right);
                have_right_plane = (I_right_Apm > 0.0);
            }

            double steer_angle_deg = std::numeric_limits<double>::quiet_NaN();
            double y_mean_pred_400mm_m = std::numeric_limits<double>::quiet_NaN();
            double y_mean_pred_500mm_m = std::numeric_limits<double>::quiet_NaN();
            double y_mean_pred_600mm_m = std::numeric_limits<double>::quiet_NaN();

            if (have_right_plane && I_ag_out_Apm > 0.0) {
                const double dx = x_right_plane - x_ag_plane;
                if (dx > 0.0) {
                    const double slope = (ymean_right - ymean_ag) / dx;
                    steer_angle_deg = rad2deg(std::atan(slope));
                    y_mean_pred_400mm_m = ymean_ag + slope * 0.4;
                    y_mean_pred_500mm_m = ymean_ag + slope * 0.5;
                    y_mean_pred_600mm_m = ymean_ag + slope * 0.6;
                }
            }

            // ---- 2.2 Collimation scan between AG exit and sample ----

            const int    Nplanes    = envi("BEAM_COL_NPLANES", 32);
            const double margin     = envd("BEAM_COL_MARGIN_M", 2.0*h);   // stay a bit away from sample wall
            const double x_start    = xa1;                // start at accel exit
            const double x_end      = g_tube_x1 - (endplate_enabled ? g_end_w : 0.0) - margin; // end just before endplate
            const double tube_inner = g_ybox - g_tube_w;  // half-height of free tube aperture

            double y_rms_max    = 0.0;
            double y_rms_c_max  = 0.0;
            double y_absmax_max = 0.0;
            bool   lost_sidewalls = false;

            for (int ip = 0; ip < Nplanes; ++ip) {
                const double s  = (Nplanes > 1) ? double(ip) / double(Nplanes - 1) : 0.0;
                const double xp = x_start + s * (x_end - x_start);

                double I_plane = 0.0, ymean = 0.0, yr = 0.0, yr_c = 0.0, ymax = 0.0;
                current_and_width_at_plane(xp, I_plane, ymean, yr, yr_c, ymax);

                if (I_plane <= 0.0)
                    continue;  // no trajectories crossing here

                if (yr > y_rms_max)      y_rms_max    = yr;
                if (yr_c > y_rms_c_max)  y_rms_c_max  = yr_c;
                if (ymax > y_absmax_max) y_absmax_max = ymax;

                // If any significant current reaches near the tube inner radius, flag it.
                if (ymax >= 0.99 * tube_inner) {
                    lost_sidewalls = true;
                }
            }

            // -------------------------------------------------------------------------
            // Moment/envelope prediction to the *physical* drift end (x_right_phys).
            // This lets us truncate the simulated drift (X_RIGHT_M) while still estimating
            // whether the beam will hit the tube wall further downstream.
            //
            // Model (simple, robust): r'' = K_eff / r, where r is the *RMS* radius and
            // K_eff is inferred from 3 RMS samples just after the accel grid exit.
            // -------------------------------------------------------------------------
            bool   env_ok        = false;
            std::string env_reason;
            double footprint_pred_m = std::numeric_limits<double>::quiet_NaN(); // predicted absmax at x_phys_end
            double clearance_pred_m = std::numeric_limits<double>::quiet_NaN();
            double x_hit_pred_m     = std::numeric_limits<double>::quiet_NaN();
            double max_abs_pred_m   = std::numeric_limits<double>::quiet_NaN();
            int    hits_tube_pred   = 0;
            double div_angle_deg    = std::numeric_limits<double>::quiet_NaN();

            double footprint_pred_400mm_m = std::numeric_limits<double>::quiet_NaN();
            double footprint_pred_500mm_m = std::numeric_limits<double>::quiet_NaN();
            double footprint_pred_600mm_m = std::numeric_limits<double>::quiet_NaN();

            const double x_phys_end = xmax_phys - g_end_w - 0.5*h; // approx "sample plane" of the full-length geometry
            const double x_targets[3] = { x_ag_plane + 0.4, x_ag_plane + 0.5, x_ag_plane + 0.6 };
            double r_targets[3] = {
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()
            };
            bool got_targets[3] = { false, false, false };
            double r_end = std::numeric_limits<double>::quiet_NaN();
            bool got_phys_end = false;

            // Minimal JSON string escaper for diagnostics.
            auto json_q = []( const std::string &s ) -> std::string {
                std::string o; o.reserve(s.size()+8);
                o.push_back('"');
                for( char c : s ) {
                    switch(c) {
                    case '\\': o += "\\\\"; break;
                    case '"':  o += "\\\""; break;
                    case '\n': o += "\\n"; break;
                    case '\r': o += "\\r"; break;
                    case '\t': o += "\\t"; break;
                    default:   o.push_back(c); break;
                    }
                }
                o.push_back('"');
                return o;
            };

            if( x_phys_end <= x_ag_plane + 2.0*h ) {
                env_reason = "x_phys_end too close to AG exit";
            } else {
                // Need a little simulated drift after AG exit to estimate curvature.
                double dx_probe = envd("ENV_PROBE_DX_M", 1.0e-3);
                dx_probe = std::max(dx_probe, 2.0*h);

                const double max_dx = (xmax - x_ag_plane - 6.0*h)/3.0;
                if( max_dx <= 2.0*h ) {
                    env_reason = "not enough simulated drift to estimate curvature";
                } else {
                    dx_probe = std::min(dx_probe, max_dx);

                    const double x0 = x_ag_plane + dx_probe;
                    const double x1 = x0 + dx_probe;
                    const double x2 = x1 + dx_probe;

                    double I0_Apm=0.0, y0_mean=0.0, r0=0.0, r0c=0.0, a0=0.0;
                    double I1_Apm=0.0, y1_mean=0.0, r1=0.0, r1c=0.0, a1=0.0;
                    double I2_Apm=0.0, y2_mean=0.0, r2=0.0, r2c=0.0, a2=0.0;

                    current_and_width_at_plane(x0, I0_Apm, y0_mean, r0, r0c, a0);
                    current_and_width_at_plane(x1, I1_Apm, y1_mean, r1, r1c, a1);
                    current_and_width_at_plane(x2, I2_Apm, y2_mean, r2, r2c, a2);

                    if( I0_Apm <= 0.0 || r0 <= 0.0 ||
                        !std::isfinite(r0) || !std::isfinite(r1) || !std::isfinite(r2) ) {
                        env_reason = "insufficient beam at probe planes";
                    } else {
                        const double rp0  = (r1 - r0) / dx_probe;
                        const double rpp0 = (r2 - 2.0*r1 + r0) / (dx_probe*dx_probe);

                        // Average divergence angle from *centered* RMS growth rate.
                        const double rp_avg = (r2c - r0c) / (2.0*dx_probe);
                        if (std::isfinite(rp_avg)) {
                            const double rp_use = std::max(0.0, rp_avg);
                            div_angle_deg = rad2deg(std::atan(rp_use));
                        }

                        double K_eff = r0 * rpp0;
                        if( !std::isfinite(K_eff) ) K_eff = 0.0;
                        if( K_eff < 0.0 ) K_eff = 0.0; // clamp: we only model divergence here

                        // Convert RMS -> "footprint" using the measured ratio absmax/RMS at x0.
                        double abs_ratio = (a0 > 0.0 && r0 > 0.0) ? (a0 / r0) : 0.0;
                        if( !(abs_ratio > 0.0) || !std::isfinite(abs_ratio) ) abs_ratio = 3.0;

                        double dx_env = envd("ENV_STEP_M", 2.0e-3);
                        dx_env = std::max(dx_env, 1.0e-4);
                        dx_env = std::min(dx_env, 1.0e-2);

                        double x  = x0;
                        double r  = r0;
                        double rp = rp0;

                        max_abs_pred_m = abs_ratio * r;
                        hits_tube_pred = 0;

                        const double x_max_target = std::max(x_phys_end,
                            std::max(x_targets[0], std::max(x_targets[1], x_targets[2])));

                        while( x < x_max_target ) {
                            const double step = std::min(dx_env, x_max_target - x);
                            const double rpp  = (r > 1.0e-12) ? (K_eff / r) : 0.0;

                            rp += rpp * step;
                            r  += rp  * step;
                            x  += step;

                            if( !std::isfinite(r) || r < 0.0 ) {
                                env_reason = "envelope diverged numerically";
                                break;
                            }

                            if (x <= x_phys_end + 1.0e-12) {
                                const double abs_now = abs_ratio * r;
                                if( abs_now > max_abs_pred_m ) max_abs_pred_m = abs_now;
                                if( !hits_tube_pred && abs_now >= tube_inner ) {
                                    hits_tube_pred = 1;
                                    x_hit_pred_m   = x;
                                }
                            }
                            if (!got_phys_end && x >= x_phys_end) {
                                got_phys_end = true;
                                r_end = r;
                            }
                            for (int ti = 0; ti < 3; ++ti) {
                                if (!got_targets[ti] && x >= x_targets[ti]) {
                                    got_targets[ti] = true;
                                    r_targets[ti] = r;
                                }
                            }
                        }

                        if( env_reason.empty() ) {
                            if (got_phys_end && std::isfinite(r_end)) {
                                footprint_pred_m = abs_ratio * r_end;
                                clearance_pred_m = tube_inner - footprint_pred_m;
                                env_ok = std::isfinite(footprint_pred_m);
                            } else {
                                env_reason = "missing envelope at physical end";
                            }
                            if (got_targets[0] && std::isfinite(r_targets[0]))
                                footprint_pred_400mm_m = abs_ratio * r_targets[0];
                            if (got_targets[1] && std::isfinite(r_targets[1]))
                                footprint_pred_500mm_m = abs_ratio * r_targets[1];
                            if (got_targets[2] && std::isfinite(r_targets[2]))
                                footprint_pred_600mm_m = abs_ratio * r_targets[2];
                        }
                    }
                }
            }

            // Write standalone envelope_pred.json (handy for debugging / plotting).
            {
                std::ofstream ej(outdir + "/envelope_pred.json");
                ej << "{\n";
                ej << "  \"ok\": " << (env_ok ? "true" : "false") << ",\n";
                ej << "  \"reason\": " << json_q(env_reason) << ",\n";
                ej << "  \"truncated_drift\": " << (truncated_drift ? "true" : "false") << ",\n";
                ej << "  \"xmax_m\": " << std::setprecision(10) << xmax << ",\n";
                ej << "  \"xmax_phys_m\": " << std::setprecision(10) << xmax_phys << ",\n";
                ej << "  \"x_ag_plane_m\": " << std::setprecision(10) << x_ag_plane << ",\n";
                ej << "  \"x_phys_end_m\": " << std::setprecision(10) << x_phys_end << ",\n";
                ej << "  \"tube_inner_m\": " << std::setprecision(10) << tube_inner << ",\n";
                ej << "  \"hits_tube_pred\": ";
                if( env_ok ) ej << (hits_tube_pred ? "true" : "false"); else ej << "null";
                ej << ",\n";
                ej << "  \"x_hit_pred_m\": ";
                if( env_ok && std::isfinite(x_hit_pred_m) ) ej << std::setprecision(10) << x_hit_pred_m; else ej << "null";
                ej << ",\n";
                ej << "  \"footprint_pred_m\": ";
                if( env_ok && std::isfinite(footprint_pred_m) ) ej << std::setprecision(10) << footprint_pred_m; else ej << "null";
                ej << ",\n";
                ej << "  \"footprint_pred_400mm_m\": ";
                if( env_ok && std::isfinite(footprint_pred_400mm_m) ) ej << std::setprecision(10) << footprint_pred_400mm_m; else ej << "null";
                ej << ",\n";
                ej << "  \"footprint_pred_500mm_m\": ";
                if( env_ok && std::isfinite(footprint_pred_500mm_m) ) ej << std::setprecision(10) << footprint_pred_500mm_m; else ej << "null";
                ej << ",\n";
                ej << "  \"footprint_pred_600mm_m\": ";
                if( env_ok && std::isfinite(footprint_pred_600mm_m) ) ej << std::setprecision(10) << footprint_pred_600mm_m; else ej << "null";
                ej << ",\n";
                ej << "  \"clearance_pred_m\": ";
                if( env_ok && std::isfinite(clearance_pred_m) ) ej << std::setprecision(10) << clearance_pred_m; else ej << "null";
                ej << ",\n";
                ej << "  \"divergence_angle_deg\": ";
                if( std::isfinite(div_angle_deg) ) ej << std::setprecision(10) << div_angle_deg; else ej << "null";
                ej << ",\n";
                ej << "  \"max_abs_pred_m\": ";
                if( env_ok && std::isfinite(max_abs_pred_m) ) ej << std::setprecision(10) << max_abs_pred_m; else ej << "null";
                ej << "\n";
                ej << "}\n";
            }
            // ---- 2.3 High-level beam summary for Stage-2 optimizer ----
            {
                // Does any beam actually reach the sample plane?
                const bool has_sample_beam = endplate_enabled ? (I_sm_Apm > 0.0) : (I_ag_out_Apm > 0.0);

                // “Good” single beam: reaches sample, does not graze tube inner wall.
                const bool good_single_beam =
                    has_sample_beam &&
                    !lost_sidewalls &&
                    (y_absmax_max < tube_inner) &&
                    (!endplate_enabled || (ymax_sm < tube_inner));

                std::ofstream mj(outdir + "/beam_metrics.json");
                mj << "{\n";
                mj << "  \"truncated_drift\": " << (truncated_drift ? "true" : "false") << ",\n";
                mj << "  \"envelope_ok\": " << (env_ok ? "true" : "false") << ",\n";
                mj << "  \"xmax_phys_m\": " << std::setprecision(10) << xmax_phys << ",\n";
                mj << "  \"x_phys_end_m\": " << std::setprecision(10) << x_phys_end << ",\n";
                mj << "  \"FOOTPRINT_PRED_M\": ";
                if( env_ok && std::isfinite(footprint_pred_m) ) mj << std::setprecision(10) << footprint_pred_m; else mj << "null";
                mj << ",\n";
                mj << "  \"FOOTPRINT_PRED_400MM_M\": ";
                if( env_ok && std::isfinite(footprint_pred_400mm_m) ) mj << std::setprecision(10) << footprint_pred_400mm_m; else mj << "null";
                mj << ",\n";
                mj << "  \"FOOTPRINT_PRED_500MM_M\": ";
                if( env_ok && std::isfinite(footprint_pred_500mm_m) ) mj << std::setprecision(10) << footprint_pred_500mm_m; else mj << "null";
                mj << ",\n";
                mj << "  \"FOOTPRINT_PRED_600MM_M\": ";
                if( env_ok && std::isfinite(footprint_pred_600mm_m) ) mj << std::setprecision(10) << footprint_pred_600mm_m; else mj << "null";
                mj << ",\n";
                mj << "  \"CLEARANCE_PRED_M\": ";
                if( env_ok && std::isfinite(clearance_pred_m) ) mj << std::setprecision(10) << clearance_pred_m; else mj << "null";
                mj << ",\n";
                mj << "  \"HITS_TUBE_PRED\": ";
                if( env_ok ) mj << (hits_tube_pred ? "true" : "false"); else mj << "null";
                mj << ",\n";
                mj << "  \"X_HIT_PRED_M\": ";
                if( env_ok && std::isfinite(x_hit_pred_m) ) mj << std::setprecision(10) << x_hit_pred_m; else mj << "null";
                mj << ",\n";
                mj << "  \"DIVERGENCE_ANGLE_DEG\": ";
                if( std::isfinite(div_angle_deg) ) mj << std::setprecision(10) << div_angle_deg; else mj << "null";
                mj << ",\n";
                mj << "  \"sample\": {\n";
                mj << "    \"I_Apm\": "      << std::setprecision(10) << I_sm_Apm  << ",\n";
                mj << "    \"I_A\": "        << std::setprecision(10) << I_sm_A    << ",\n";
                mj << "    \"y_rms_m\": "    << std::setprecision(10) << yr_sm     << ",\n";
                mj << "    \"y_absmax_m\": " << std::setprecision(10) << ymax_sm   << "\n";
                mj << "  },\n";
                mj << "  \"currents\": {\n";
                mj << "    \"I_pg_in_Apm\": "   << std::setprecision(10) << I_pg_in_Apm  << ",\n";
                mj << "    \"I_ag_out_Apm\": "  << std::setprecision(10) << I_ag_out_Apm << ",\n";
                mj << "    \"I_pg_in_A\": "     << std::setprecision(10) << I_pg_in_A    << ",\n";
                mj << "    \"I_ag_out_A\": "    << std::setprecision(10) << I_ag_out_A   << "\n";
                mj << "  },\n";
                mj << "  \"collimation\": {\n";
                mj << "    \"y_rms_max_m\": "      << std::setprecision(10) << y_rms_max    << ",\n";
                mj << "    \"y_rms_c_max_m\": "    << std::setprecision(10) << y_rms_c_max  << ",\n";
                mj << "    \"y_absmax_max_m\": "   << std::setprecision(10) << y_absmax_max << ",\n";
                mj << "    \"tube_inner_half_m\": " << std::setprecision(10) << tube_inner   << ",\n";
                mj << "    \"lost_to_sidewalls\": " << (lost_sidewalls   ? "true" : "false") << ",\n";
                mj << "    \"good_single_beam\": "  << (good_single_beam ? "true" : "false") << ",\n";
                mj << "    \"has_sample_beam\": "   << (has_sample_beam  ? "true" : "false") << "\n";
                mj << "  }\n";
                mj << ",\n";
                mj << "  \"deflection\": {\n";
                mj << "    \"x_ag_plane_m\": " << std::setprecision(10) << x_ag_plane << ",\n";
                mj << "    \"x_right_plane_m\": " << std::setprecision(10) << x_right_plane << ",\n";
                mj << "    \"y_mean_ag_m\": " << std::setprecision(10) << ymean_ag << ",\n";
                mj << "    \"y_mean_right_m\": " << std::setprecision(10) << ymean_right << ",\n";
                mj << "    \"y_rms_c_ag_m\": " << std::setprecision(10) << yr_ag_c << ",\n";
                mj << "    \"y_rms_c_right_m\": " << std::setprecision(10) << yr_right_c << ",\n";
                mj << "    \"steer_angle_deg\": ";
                if (std::isfinite(steer_angle_deg)) mj << std::setprecision(10) << steer_angle_deg; else mj << "null";
                mj << ",\n";
                mj << "    \"y_mean_pred_400mm_m\": ";
                if (std::isfinite(y_mean_pred_400mm_m)) mj << std::setprecision(10) << y_mean_pred_400mm_m; else mj << "null";
                mj << ",\n";
                mj << "    \"y_mean_pred_500mm_m\": ";
                if (std::isfinite(y_mean_pred_500mm_m)) mj << std::setprecision(10) << y_mean_pred_500mm_m; else mj << "null";
                mj << ",\n";
                mj << "    \"y_mean_pred_600mm_m\": ";
                if (std::isfinite(y_mean_pred_600mm_m)) mj << std::setprecision(10) << y_mean_pred_600mm_m; else mj << "null";
                mj << "\n";
                mj << "  }\n";
                // Perveance numbers (geometric from simulated current, plus CL reference and normalized perveance)
//
// IMPORTANT NOTE (truncated drift):
// If you shrink X_RIGHT_M to drop most of the downstream drift tube, the "sample plane" diagnostic
// may land outside the physically modeled region (or very near an open/Neumann boundary).
// In that case it is common to see I_sm_A -> 0 even when the extracted beam at the accelerator exit
// is healthy. For perveance tracking we therefore *default* to AG-exit current and also report the
// "true" sample-plane perveance separately when available.
const double P_geom_pg_A_per_V32 = (V32 > 0.0) ? (I_pg_in_A  / V32) : 0.0;
const double P_geom_ag_A_per_V32 = (V32 > 0.0) ? (I_ag_out_A / V32) : 0.0;

// True sample-plane perveance (only meaningful when the sample plane is physically present)
const double P_geom_sm_true_A_per_V32 = (V32 > 0.0) ? (I_sm_A / V32) : 0.0;

// Back-compat: the "sm" perveance fields now follow AG-exit current to avoid spurious zeros
const double P_geom_sm_A_per_V32 = P_geom_ag_A_per_V32;

const double Pnorm_pg = (P_CL_A_per_V32 > 0.0) ? (P_geom_pg_A_per_V32 / P_CL_A_per_V32) : 0.0;
const double Pnorm_ag = (P_CL_A_per_V32 > 0.0) ? (P_geom_ag_A_per_V32 / P_CL_A_per_V32) : 0.0;
const double Pnorm_sm_true = (P_CL_A_per_V32 > 0.0) ? (P_geom_sm_true_A_per_V32 / P_CL_A_per_V32) : 0.0;
const double Pnorm_sm = Pnorm_ag;

const double I_sys_pg_A = I_pg_in_A  * N_beamlets_est;
const double I_sys_ag_A = I_ag_out_A * N_beamlets_est;
const double I_sys_sm_true_A = I_sm_A * N_beamlets_est;

// Back-compat: I_sample_A reported in the perveance block follows AG-exit current
const double I_sys_sm_A = I_sys_ag_A;

const double P_sys_geom_pg_A_per_V32 = P_geom_pg_A_per_V32 * N_beamlets_est;
const double P_sys_geom_ag_A_per_V32 = P_geom_ag_A_per_V32 * N_beamlets_est;
const double P_sys_geom_sm_true_A_per_V32 = P_geom_sm_true_A_per_V32 * N_beamlets_est;

// Back-compat: system "sm" perveance follows AG-exit current
const double P_sys_geom_sm_A_per_V32 = P_sys_geom_ag_A_per_V32;

const double P_sys_norm_pg = (P_CL_sys_A_per_V32 > 0.0) ? (P_sys_geom_pg_A_per_V32 / P_CL_sys_A_per_V32) : 0.0;
const double P_sys_norm_ag = (P_CL_sys_A_per_V32 > 0.0) ? (P_sys_geom_ag_A_per_V32 / P_CL_sys_A_per_V32) : 0.0;
const double P_sys_norm_sm_true = (P_CL_sys_A_per_V32 > 0.0) ? (P_sys_geom_sm_true_A_per_V32 / P_CL_sys_A_per_V32) : 0.0;
const double P_sys_norm_sm = P_sys_norm_ag;


                mj << ",\n";
                mj << "  \"perveance\": {\n";
                mj << "    \"V_accel_V\": " << std::setprecision(10) << V_accel_V << ",\n";
                mj << "    \"V32\": "       << std::setprecision(10) << V32       << ",\n";
                mj << "    \"d_eff_m\": "   << std::setprecision(10) << d_eff_m   << ",\n";
                mj << "    \"P_CL_A_per_V32\": " << std::setprecision(10) << P_CL_A_per_V32 << ",\n";
                mj << "    \"I_CL_A\": "         << std::setprecision(10) << I_CL_A         << ",\n";
                mj << "    \"N_beamlets_est\": " << std::setprecision(10) << N_beamlets_est << ",\n";
                mj << "    \"P_CL_sys_A_per_V32\": " << std::setprecision(10) << P_CL_sys_A_per_V32 << ",\n";
                mj << "    \"I_CL_sys_A\": "         << std::setprecision(10) << I_CL_sys_A         << ",\n";

                mj << "    \"beamlet\": {\n";
                mj << "      \"P_geom_pg_A_per_V32\": " << std::setprecision(10) << P_geom_pg_A_per_V32 << ",\n";
                mj << "      \"P_geom_ag_A_per_V32\": " << std::setprecision(10) << P_geom_ag_A_per_V32 << ",\n";
                mj << "      \"P_geom_sm_A_per_V32\": " << std::setprecision(10) << P_geom_sm_A_per_V32 << ",\n";
                mj << "      \"P_geom_sm_true_A_per_V32\": " << std::setprecision(10) << P_geom_sm_true_A_per_V32 << ",\n";
                mj << "      \"P_norm_pg\": " << std::setprecision(10) << Pnorm_pg << ",\n";
                mj << "      \"P_norm_ag\": " << std::setprecision(10) << Pnorm_ag << ",\n";
                mj << "      \"P_norm_sm\": " << std::setprecision(10) << Pnorm_sm << ",\n";
                mj << "      \"P_norm_sm_true\": " << std::setprecision(10) << Pnorm_sm_true << "\n";
                mj << "    },\n";

                mj << "    \"system\": {\n";
                mj << "      \"I_pg_in_A\": "   << std::setprecision(10) << I_sys_pg_A << ",\n";
                mj << "      \"I_ag_out_A\": "  << std::setprecision(10) << I_sys_ag_A << ",\n";
                mj << "      \"I_sample_A\": "  << std::setprecision(10) << I_sys_sm_A << ",\n";
                mj << "      \"I_sample_true_A\": "  << std::setprecision(10) << I_sys_sm_true_A << ",\n";
                mj << "      \"P_geom_pg_A_per_V32\": " << std::setprecision(10) << P_sys_geom_pg_A_per_V32 << ",\n";
                mj << "      \"P_geom_ag_A_per_V32\": " << std::setprecision(10) << P_sys_geom_ag_A_per_V32 << ",\n";
                mj << "      \"P_geom_sm_A_per_V32\": " << std::setprecision(10) << P_sys_geom_sm_A_per_V32 << ",\n";
                mj << "      \"P_geom_sm_true_A_per_V32\": " << std::setprecision(10) << P_sys_geom_sm_true_A_per_V32 << ",\n";
                mj << "      \"P_norm_pg\": " << std::setprecision(10) << P_sys_norm_pg << ",\n";
                mj << "      \"P_norm_ag\": " << std::setprecision(10) << P_sys_norm_ag << ",\n";
                mj << "      \"P_norm_sm\": " << std::setprecision(10) << P_sys_norm_sm << ",\n";
                mj << "      \"P_norm_sm_true\": " << std::setprecision(10) << P_sys_norm_sm_true << "\n";
                mj << "    }\n";

                mj << "  }\n";
                mj << "}\n";
            }


            if (lost_sidewalls) {
                std::fprintf(stderr,
                            "[beam] ERROR: beam reaches tube inner wall (|y| >= %.3e m) "
                            "between AG exit and sample; significant loss to sidewalls likely.\n",
                            tube_inner);
            }

        }

        t_after_solve = SteadyClock::now();

        t_after_diag = SteadyClock::now();

        if (writepng) {
            // Full view
	    const int FONT_PX = envi("PLOT_FONT_PX", 22);
            GeomPlotter gp(geom);
            gp.set_size(1600, 800);
            gp.set_view(VIEW_XY, -1);
            gp.set_epot(&epot);
	    gp.set_font_size(FONT_PX);
            gp.set_efield(&efield);
            if (use_pdb) {
                gp.set_particle_database(&pdb);  // overlay ion trajectories when ENABLE_IONS=1
            }
            gp.plot_png(png.c_str());
            std::printf("wrote PNG (with epot/efield%s): %s\n",
                        use_pdb ? " + beam" : "", png.c_str());

            // Close-up around the grids
            const double pad = envd("CLOSE_PAD_M", 5e-4);  // 0.5 mm default
            double x0 = xs0 - pad;
            double x1 = xa1 + pad;

            // account for offset + chamfer widening (centered between screen/accel offsets)
            double max_delta2 = std::max({ g_scr_du*g_scr_mu, g_scr_dd*g_scr_md,
                                           g_acc_du*g_acc_mu, g_acc_dd*g_acc_md });
            const double y_center = 0.5 * (scr_off_y + off_y);
            const double y_span = std::max(std::fabs(scr_off_y - y_center),
                                           std::fabs(off_y - y_center));
            double ymax_local = y_span + a + max_delta2 + pad;
            double y0 = y_center - ymax_local, y1 = y_center + ymax_local;

            // aspect-preserving size: 800 px tall, width scaled by xr/yr
            const int   BASE = envi("CLOSE_BASE_PX", 800);
            const double xr = x1 - x0, yr = y1 - y0;
            const int   H = BASE;
            const int   W = std::max(1, int(std::lround(BASE * (xr/yr))));

            GeomPlotter gpc(geom);
            gpc.set_view(VIEW_XY, -1);
            gpc.set_epot(&epot);
            gpc.set_efield(&efield);
            if (use_pdb) {
                gpc.set_particle_database(&pdb);
            }
            gpc.set_ranges(x0, y0, x1, y1);
	    gpc.set_font_size(FONT_PX);
            gpc.set_size(W, H);

            // write close-up next to main PNG in results dir
            std::string close = add_suffix_before_ext(png, "_closeup");
            gpc.plot_png(close.c_str());
            std::printf("wrote PNG (close-up%s): %s\n",
                        use_pdb ? " + beam" : "", close.c_str());
        }

        t_after_png = SteadyClock::now();

    } else if (writepng) {
    	const int FONT_PX = envi("PLOT_FONT_PX", 22);

        // Geometry-only PNGs (no field solve)
        GeomPlotter gp(geom);
        gp.set_size(1600, 800);
        gp.set_view(VIEW_XY, -1);
	gp.set_font_size(FONT_PX);
        gp.plot_png(png.c_str());
        std::printf("wrote PNG: %s\n", png.c_str());

        // Close-up
        const double pad = envd("CLOSE_PAD_M", 5e-4);
        double x0 = xs0 - pad;
        double x1 = xa1 + pad;
        double max_delta = std::max({ g_scr_du*g_scr_mu, g_scr_dd*g_scr_md,
                                      g_acc_du*g_acc_mu, g_acc_dd*g_acc_md });
        const double y_center = 0.5 * (scr_off_y + off_y);
        const double y_span = std::max(std::fabs(scr_off_y - y_center),
                                       std::fabs(off_y - y_center));
        double ymax_local = y_span + a + max_delta + pad;
        double y0 = y_center - ymax_local, y1 = y_center + ymax_local;

        const int   BASE = envi("CLOSE_BASE_PX", 800);
        const double xr = x1 - x0, yr = y1 - y0;
        const int   H = BASE;
        const int   W = std::max(1, int(std::lround(BASE * (xr/yr))));

        GeomPlotter gpc(geom);
        gpc.set_view(VIEW_XY, -1);
        gpc.set_ranges(x0, y0, x1, y1);
        gpc.set_size(W, H);
	gpc.set_font_size(FONT_PX);

        // write close-up next to main PNG in results dir
        std::string close = add_suffix_before_ext(png, "_closeup");
        gpc.plot_png(close.c_str());
        std::printf("wrote PNG (close-up): %s\n", close.c_str());
        t_after_png = SteadyClock::now();
    }




// --- write timing.json (safe string literals) ---
try {
    t_end = SteadyClock::now();
    std::ofstream tj(outdir + "/timing.json");
    tj << "{\n";
    tj << "  \"wall_total_s\": "  << std::setprecision(10) << sec_since(t_start, t_end) << ",\n";
    tj << "  \"mesh_build_s\": "  << std::setprecision(10) << sec_since(t_start, t_after_mesh) << ",\n";
    tj << "  \"solve_block_s\": " << std::setprecision(10) << sec_since(t_after_mesh, t_after_solve) << ",\n";
    tj << "  \"diagnostics_s\": " << std::setprecision(10) << sec_since(t_after_solve, t_after_diag) << ",\n";
    tj << "  \"png_s\": "        << std::setprecision(10) << sec_since(t_after_diag, t_after_png) << "\n";
    tj << "}\n";
} catch (...) {
    // don't fail the run if timing write fails
}


    return 0;
}

catch (Error &e) {
    e.print_error_message(std::cerr, true);
    return 1;
}
