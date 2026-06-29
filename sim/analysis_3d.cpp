// ============================================================================
//  analysis_3d.cpp
//  Interactive 3-D viewer for a finished grid_3d_curved run, using IBSimu's
//  GTKPlotter. Loads the saved state and opens a rotatable geometry+trajectory
//  window (isometric views of the beamline).
//
//  Requires a display: run on a workstation with IBSimu+GTK, or over `ssh -X`
//  from the cluster (works, just laggy when rotating).
//
//  Usage:
//      ./analysis_3d results_3d/geom.dat results_3d/epot.dat results_3d/pdb.dat
//  Build:
//      make -f Makefile.3d analysis_3d
// ============================================================================
#include <fstream>
#include <iostream>
#include "geometry.hpp"
#include "epot_efield.hpp"
#include "meshscalarfield.hpp"
#include "particledatabase.hpp"
#include "gtkplotter.hpp"
#include "error.hpp"
#include "ibsimu.hpp"

void simu( int argc, char **argv )
{
    if( argc < 4 ) {
        std::cerr << "Usage: analysis_3d geom.dat epot.dat pdb.dat\n";
        throw( Error( ERROR_LOCATION, "missing input files" ) );
    }

    std::ifstream is_geom( argv[1] );
    if( !is_geom.good() ) throw( Error( ERROR_LOCATION, (std::string)"cannot open "+argv[1] ) );
    Geometry geom( is_geom );
    is_geom.close();
    geom.build_surface();

    std::ifstream is_epot( argv[2] );
    if( !is_epot.good() ) throw( Error( ERROR_LOCATION, (std::string)"cannot open "+argv[2] ) );
    EpotField epot( is_epot, geom );
    is_epot.close();

    EpotEfield efield( epot );
    field_extrpl_e ex[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                             FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                             FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( ex );

    std::ifstream is_pdb( argv[3] );
    if( !is_pdb.good() ) throw( Error( ERROR_LOCATION, (std::string)"cannot open "+argv[3] ) );
    ParticleDataBase3D pdb( is_pdb, geom );
    is_pdb.close();

    // trajectory density (nice for a 3-D beam-envelope overlay)
    MeshScalarField tdens( geom );
    pdb.build_trajectory_density_field( tdens );

    GTKPlotter plotter( &argc, &argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.set_efield( &efield );
    plotter.set_trajdens( &tdens );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();
}

int main( int argc, char **argv )
{
    try {
        ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
        ibsimu.set_thread_count( 1 );
        simu( argc, argv );
    } catch( Error e ) {
        e.print_error_message( ibsimu.message(0) );
        return 1;
    }
    return 0;
}
