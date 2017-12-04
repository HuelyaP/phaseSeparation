/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Peter Weisbrod
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

/* phaseSeparation2d.cpp:
 * In this example the simulation is initialized with a given
 * density plus a small random number all over the domain. This
 * condition is unstable and leads to liquid-vapor phase separation.
 * Boundaries are assumed to be periodic. This example shows the
 * usage of multiphase flow.
 */


#include "olb2D.h"
#include "olb2D.hh"   // use only generic version!
#include <cstdlib>
#include <iostream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;
#define DESCRIPTOR ShanChenDynOmegaForcedD2Q9Descriptor


// Parameters for the simulation setup
const int maxIter  = 5000;
const int nx   = 100;
const int ny   = 100;


// Stores geometry information in form of material numbers
void prepareGeometry(LBconverter<T> const& converter, 
							SuperGeometry2D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  // Sets material number for fluid
  superGeometry.rename( 0,1 );

  // Vectors for Geometry Settings
  Vector<T, 2> origin1(-.5, -.5);
  Vector<T, 2> extend1(.5 + nx, .5+ny/2);

  //Indicator functors
  IndicatorCuboid2D<T> bottom(extend1, origin1);
  superGeometry.rename(1, 2, bottom);

  // Removes all not needed boundary voxels outside the surface
  //superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  //superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(LBconverter<T> const& converter,
							SuperLattice2D<T, DESCRIPTOR>& sLattice,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics1,
                     SuperGeometry2D<T>& superGeometry ) {
	//Dynamics 0
	sLattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>());

  //BulkDynamics
  sLattice.defineDynamics( superGeometry, 1, &bulkDynamics1 );
  sLattice.defineDynamics(superGeometry, 2, &bulkDynamics1);

  // Initial conditions
  //AnalyticalConst2D<T,T> noise( 2. );
  std::vector<T> v( 2,T() );
  AnalyticalConst2D<T,T> zeroVelocity( v );
  //AnalyticalConst2D<T,T> rhoUpper( 2.659 ); //werte aus timm krüger S375
  AnalyticalConst2D<T, T> rhoUpper(8.);
  //AnalyticalConst2D<T, T> rhoBottom(.056);
  AnalyticalConst2D<T, T> rhoBottom(.1);
  AnalyticalRandom2D<T,T> random;
  AnalyticalIdentity2D<T,T> upper(rhoUpper );
  AnalyticalIdentity2D<T, T> bottom(rhoBottom);

  //just a test

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( superGeometry, 1, upper, zeroVelocity );
  sLattice.iniEquilibrium( superGeometry, 1, upper, zeroVelocity );
  sLattice.defineRhoU(superGeometry, 2, bottom, zeroVelocity);
  sLattice.iniEquilibrium(superGeometry, 2, bottom, zeroVelocity);

  // Make the lattice ready for simulation
  sLattice.initialize();
}

void saveRhoEq(double fluid, double vapor) {
	//Create a file object
	ofstream rhoEq;
	//Open the file
	rhoEq.open("rhoEq.txt", ios_base::out | ios_base::app);
	//Error handling
	if (!rhoEq) {
		cerr << "Fehler im File eqRho" << endl;
		getchar();
	}
	rhoEq << fluid << "," << vapor << "\n" << endl;
	//Close the file
	rhoEq.close();
}

// Output to console and files
void getResults( SuperLattice2D<T, DESCRIPTOR>& sLattice, 
						LBconverter<T> const& converter, int iT,
						SuperGeometry2D<T>& superGeometry, Timer<T>& timer) {

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter2D<T> vtmWriter( "phaseSeparation2d" );
  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(sLattice, converter);
  SuperLatticeDensity2D<T, DESCRIPTOR> density( sLattice );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( density );

  const int vtkIter  = 20;
  const int statIter = 20;

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
  }

  // Writes the vtk files
  if ( iT%vtkIter==0 ) {
    clout << "Writing VTK..." << std::endl;
    vtmWriter.write( iT );

    BlockLatticeReduction2D<T, DESCRIPTOR> planeReduction( density );
    BlockGifWriter<T> gifWriter;
    gifWriter.write( planeReduction, iT, "density" );
  }

  // Writes output on the console
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT, converter.physTime(iT));
  }

  /***********************************/
  // Output for x-velocity along y-position at the last time step
  if (iT == maxIter - 1) {
	  // Gives access to density information on lattice in lattice units
	  SuperLatticeDensity2D<T, DESCRIPTOR> densityField(sLattice);
	  // Interpolation functor with densityField information - Ist das nötig ???
	  AnalyticalFfromSuperLatticeF2D<T, DESCRIPTOR> interpolation(densityField, true);

	  //Vector for density_simulation
	  Vector<T, nx> density_simulation;
	  // Gnuplot interface to create plots
	  static Gnuplot<T> gplot("interfaceTest_");
	  for (int nY = 0; nY < nx; ++nY) {

		  T position[2] = { 50, nY / 1. }; 
		  T density[2] = { T(), T() }; 
			// Interpolate densityField at "position" and save it in "density"
		  interpolation(density, position);
		  // y-values: value of density (in x-direction) in "density_simulation" for every position "nY"
		  density_simulation[nY] = density[0];
		  //Extract density of the fluid an of the vapor
		  double vapor;
		  double fluid;
		  if (nY == 5) {
			  fluid = density_simulation[nY];
		  }
		  if (nY == nx / 2) {
			  vapor = density_simulation[nY];
		  }
		  if (nY == nx - 1) {
			  saveRhoEq(fluid, vapor);
		  }
		  // Set data for plot output
		  gplot.setData(position[1], density_simulation[nY], { "simulated" });
	  }
	  // Create PNG file
	  gplot.writePNG();
  }//end if
}

int main( int argc, char *argv[] ) {

  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  //clout.setMultiOutput(true);

	//converter----------------------------------------------
	LBconverter<T> converter(
	( int )	2,										//dim
	( T )		1./ nx,								//latticeL
	( T )		0.1,									//latticeU,
	( T )		1.16 *	pow(10, -7),			//charNu ??
	( T )		0.03,									//charL = 1,
	( T )		6.96*pow(10, -6),					//charU = 1,	
	( T )		107.									//	charRho
	);
	converter.print();
	writeLogFile(converter, "phaseSeperation");

	//const T omega1 = converter.getOmega();
	const T G      = -6;
	//const T G = -0.0649; //für T = 606,14K
	// const T omega1 = converter.getOmega();

  // === 2rd Step: Prepare Geometry ===

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidGeometry2D<T> cuboidGeometry( 0, 0, 1, nx, ny, noOfCuboids );

  // Periodic boundaries in x- and y-direction
  cuboidGeometry.setPeriodicity( true, true );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry2D<T> superGeometry( cuboidGeometry,loadBalancer,2 );

  prepareGeometry( converter, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice2D<T, DESCRIPTOR> sLattice( superGeometry );

  ForcedShanChenBGKdynamics<T, DESCRIPTOR> bulkDynamics1 (
    converter.getOmega(), instances::getExternalVelocityMomenta<T,DESCRIPTOR>() );

  std::vector<T> rho0;
  rho0.push_back( 1 );
  rho0.push_back( 1 );
  //ShanChen93<T,T> interactionPotential;
  
  PengRobinson<T, T> interactionPotential(G);
  //hier PengRobinson
  ShanChenForcedSingleComponentGenerator2D<T,DESCRIPTOR> coupling( G,rho0,interactionPotential );

  sLattice.addLatticeCoupling( superGeometry, 1, coupling, sLattice );
  sLattice.addLatticeCoupling(superGeometry, 2, coupling, sLattice);

  prepareLattice( converter, sLattice, bulkDynamics1, superGeometry );

  // === 4th Step: Main Loop ===
  //Konvergenzkriterium
  //T interval = 20;
  //T epsilon = 0.01;
  //util::ValueTracer<T> converge(interval, epsilon);

  int iT = 0;
  clout << "starting simulation..." << endl;
  Timer<T> timer(maxIter, superGeometry.getStatistics().getNvoxel());
  timer.start();
  for (iT = 0; iT < maxIter; ++iT) {

	  //Falls System konvergiert
	  //if (converge.hasConverged()) {
		 // clout << "Simulation converged." << endl;
		 // getResults(sLattice, converter, iT, superGeometry, timer, converge.hasConverged()); //hier
		 // break;
	  //}

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    // in this application no boundary conditions have to be adjusted

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    sLattice.communicate();
    sLattice.executeCoupling();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer);
	 //converge.takeValue(sLattice.getStatistics().getAverageEnergy(), true);
  }

  timer.stop();
  timer.printSummary();
}

