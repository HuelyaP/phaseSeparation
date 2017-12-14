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
#include <fstream>
#include <string>
#include <sstream>
#include <math.h> /* pow */

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;
#define DESCRIPTOR ShanChenDynOmegaForcedD2Q9Descriptor

// Parameters for the simulation setup
const int maxIter = 10000;
//const int maxIter = 1; //zur Kontrolle, ob Interface aufgebaut wurde
//Auflösung
const int N = 10000;
const double nx = 0.01;
const double ny = 0.01;
//const T radius = 0.005; //0.002

//std::string rad = boost::lexical_cast<std::string>(radius);

string toString(double x) {
	ostringstream x_convert;
	x_convert << x;
	return x_convert.str();
}



// Stores geometry information in form of material numbers
void prepareGeometry(LBconverter<T> const& converter, 
							SuperGeometry2D<T>& superGeometry, T radius) {

	OstreamManager clout(std::cout, "prepareGeometry");
	clout << "Prepare Geometry ..." << std::endl;

	Vector<T, 2> origin;
	Vector<T, 2> extend(nx, ny);
	IndicatorCuboid2D<T> volume(extend, origin);

	//T radius = 0.2 * nx;
	Vector<T, 2> center(nx/2, ny/2);
	//SmoothIndicatorCircle2D<T, T> circle(center, radius, 1, 10);
	IndicatorCircle2D<T> circle_inner(center, radius);

	// remove step from channel
	IndicatorIdentity2D<T> channelIdent(volume - circle_inner);

	// Sets material number for fluid and boundary
	superGeometry.rename(0, 1, channelIdent); //channel
	//superGeometry.rename(0, 1);
	
	// Removes all not needed boundary voxels outside the surface
	superGeometry.clean();

	// Removes all not needed boundary voxels inside the surface
	superGeometry.innerClean();
	superGeometry.checkForErrors();
	superGeometry.print();

	clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(LBconverter<T> const& converter, SuperLattice2D<T, DESCRIPTOR>& sLattice,
	Dynamics<T, DESCRIPTOR>& bulkDynamics1,
	SuperGeometry2D<T>& superGeometry) {

	OstreamManager clout(std::cout, "prepareLattice");
	clout << "Prepare Lattice ..." << std::endl;

	// Material=1,2 -->bulk dynamics, above
	//sLattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>());
	sLattice.defineDynamics(superGeometry, 0, &bulkDynamics1);
	sLattice.defineDynamics(superGeometry, 1, &bulkDynamics1);
	

	// Initial conditions
	std::vector<T> v(2, T());
	AnalyticalConst2D<T, T> zeroVelocity(v);
	//convert: physical pressure in Pa to lattice density latticeRho = physical pressure / pressureFactor * 3 - 1
	
	//T p0 = 8.*converter.getLatticeNu()*converter.getLatticeU()*Lx / (Ly*Ly);
	//AnalyticalLinear2D<T, T> rho(-p0 / lx*DESCRIPTOR<T>::invCs2, 0, p0*DESCRIPTOR<T>::invCs2 + 1);

	AnalyticalConst2D<T, T> rhoLiquid(5.9); //PengRobinson
	//AnalyticalConst2D<T, T> rhoLiquid(2.64); //SC93
	AnalyticalConst2D<T, T> rhoVapor(0.2); //PengRobinson
	//AnalyticalConst2D<T, T> rhoVapor(0.073); //SC93

	//T radius = 0.2 * nx;
	//Vector<T, 2> center(nx / 2, ny / 2);
	//SmoothIndicatorCircle2D<T, T> circle(center, radius, 1, 10);
	//AnalyticalIdentity2D<T, T> initialRho(rhoLiquid + (rhoVapor - rhoLiquid) * circle);

	//sLattice.defineRhoU(superGeometry, 1, initialRho, zeroVelocity);
	//sLattice.iniEquilibrium(superGeometry, 1, initialRho, zeroVelocity);

	// Initialize all values of distribution functions to their local equilibrium
	sLattice.defineRhoU(superGeometry, 1, rhoLiquid, zeroVelocity);
	sLattice.iniEquilibrium(superGeometry, 1, rhoLiquid, zeroVelocity);
	sLattice.defineRhoU(superGeometry, 0, rhoVapor, zeroVelocity);
	sLattice.iniEquilibrium(superGeometry, 0, rhoVapor, zeroVelocity);
	//sLattice.defineRhoU(superGeometry, 2, rhoMix, zeroVelocity);
	//sLattice.iniEquilibrium(superGeometry, 2, rhoMix, zeroVelocity);

	// Make the lattice ready for simulation
	sLattice.initialize();

	clout << "Prepare Lattice ... OK" << std::endl;
}

void surfaceTension(T radius, double inside, double outside, LBconverter<T> const& converter) {

	//Create a file object
	ofstream surfaceTension;

	//Open the file
	surfaceTension.open("surfaceTension.txt", ios_base::out | ios_base::app);
	
	//Error handling
	if (!surfaceTension) {
		cerr << "Fehler im File surfaceTension" << endl;
		getchar();
	}
	
	//Compute pressure in SI-Units
	double insideReal = converter.physPressureFromRho(inside);
	double outsideReal = converter.physPressureFromRho(outside);
	double differenceReal = insideReal - outsideReal;

	//Compute pressure in Lattice-Units
	double insideP = inside / olb::descriptors::ShanChenForcedD2Q9Descriptor< T >::invCs2;
	double outsideP = outside / olb::descriptors::ShanChenForcedD2Q9Descriptor< T >::invCs2;
	double difference = insideP - outsideP;

	//Insert the parameters 
	surfaceTension << radius << "," << differenceReal << "," << insideReal << "," << outsideReal << "," << difference << "," << insideP << "," << outsideP << "," << inside << "," << outside << "\n" << endl;
	
	//Close the file
	surfaceTension.close();
}

	// Output to console and files
	void getResults(SuperLattice2D<T, DESCRIPTOR>& sLattice, 
						LBconverter<T> const& converter, int iT,
						SuperGeometry2D<T>& superGeometry, Timer<T>& timer,
						bool converged, T radius){

	OstreamManager clout(std::cout, "getResults");
	string radiusStr = toString(radius);

	SuperVTMwriter2D<T> vtmWriter("phaseSeparation2d" + radiusStr);
	SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(sLattice, converter );
	SuperLatticeDensity2D<T, DESCRIPTOR> density(sLattice);
	vtmWriter.addFunctor(velocity);
	vtmWriter.addFunctor(density);

	const int vtkIter = 20;
	const int statIter = 20;

	if (iT == 0) {
		// Writes the geometry, cuboid no. and rank no. as vti file for visualization
		SuperLatticeGeometry2D<T, DESCRIPTOR> geometry(sLattice, superGeometry);
		SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid(sLattice);
		SuperLatticeRank2D<T, DESCRIPTOR> rank(sLattice);
		vtmWriter.write(geometry);
		vtmWriter.write(cuboid);
		vtmWriter.write(rank);

		vtmWriter.createMasterFile();
	}

	// Writes the vtk files
	if (iT%vtkIter == 0) {
		clout << "Writing VTK..." << std::endl;
		vtmWriter.write(iT);

		BlockLatticeReduction2D<T, DESCRIPTOR> planeReduction(density);
		BlockGifWriter<T> gifWriter;
		gifWriter.write(planeReduction, iT, "density");
	}

	// Writes output on the console
	if (iT%statIter == 0) {
		// Timer console output
		timer.update(iT);
		timer.printStep();

		// Lattice statistics console output
		sLattice.getStatistics().print(iT, converter.physTime(iT));
	}
	//hier doppelt??????
	/***********************************/
	// Output for x-velocity along y-position at the last time step
	if ((iT == maxIter-1) || converged) {
		// Gives access to density information on lattice in lattice units
		SuperLatticeDensity2D<T, DESCRIPTOR> densityField(sLattice);
		// Interpolation functor with densityField information - Ist das nötig ???
		AnalyticalFfromSuperLatticeF2D<T, DESCRIPTOR> interpolation(densityField, true); 
	
		//Vector for density_simulation
		Vector<T, N> density_simulation;
		//Wozu???
		//Vector<int, 17> y_coord({ 128, 125, 124, 123, 122, 109, 94, 79, 64, 58, 36, 22, 13, 9, 8, 7, 0 });
		// Gnuplot interface to create plots
		
		static Gnuplot<T> gplot("centerDensityX_" + radiusStr);
		double inside;
		double outside;
		int step = 120;
		for (int nY = 0; nY < step; ++nY) {
			//std::cout << "Juhu! " << std::cout;
			// 17 data points evenly distributed between 0 and 1 (height)
			//x values = position[1] Gerade geht durch Punkte(x|50)
			T position[2] = { nx/2, nY*nx/step }; //???
			T density[2] = { T(), T() }; //?
			// Interpolate densityField at "position" and save it in "density"
			interpolation(density, position);
			// y-values: value of density (in x-direction) in "density_simulation" for every position "nY"
			density_simulation[nY] = density[0];
			//Extract oustide and inside density for surface tension
			//std::cout << converter.physRho(density_simulation[nY]) << std:: endl;
			
			if (nY*nx / step == 0.001) { // nY*nx/step == 0.001
				//std::cout << "oustide: " << nY <<std::endl;
				
				outside = density_simulation[nY];
			}
			if (nY*nx / step == 0.005) {
				//std::cout << "inside: " << nY <<std::endl;
				inside = density_simulation[nY];
			}
			// Set data for plot output	
			//std::cout << "Ny: " << nY << "density: " << density_simulation[nY] << std::endl;
			gplot.setData(position[1], density_simulation[nY], { "simulated"});
		}
		std::cout << "here" << std::endl;
		surfaceTension(radius, inside, outside, converter);
		// Create PNG file
		gplot.writePNG();
	}//end if

}

int main(int argc, char *argv[]) {

	
	//T radius = 0.0001;
			//check arguments
			/*if (argc != 2)
			{
				cerr << "Need XML-File" << endl;
				return 1;
			}*/

			// === 1st Step: Initialization ===
		olbInit(&argc, &argv);
		OstreamManager clout(std::cout, "main");

		int iter = 0;
		while (iter < 30) {
			T radius = 0.0001 + iter * 0.0001;
		/*string fName( argv[1] );
		XMLreader config(fName);

		std::string olbdir, outputdir;
		config["Application"]["OlbDir"].read(olbdir);
		config["Output"]["OutputDir"].read(outputdir);
		singleton::directories().setOlbDir(olbdir);
		singleton::directories().setOutputDir(outputdir);*/

		singleton::directories().setOutputDir("./tmp/");

		// display messages from every single mpi process
		//clout.setMultiOutput(true);

		LBconverter<T> converter(
			(int)2,											// dim					dimension of the domain 
			(T)	0.0001,										// latticeL_			length of a lattice cell in meter (proportional to Knudsen number) 
			(T)	0.1 * 100 / N,								// latticeU_ //0.1	velocity in dimensionless lattice units (proportional to Mach number) 
			(T)   1.16 *	pow(10, -7),					// charNu_ 				kinematic viscosity in m^2/s
			(T)   0.01,										// charL_ = 1,			characteristical length in meter
			(T)	6.96*pow(10, -6),						// charU					characteristical speed in m/s
			(T)	107.										//	charRho				density factor in kg/m^d (latticeRho can be multplied by this factor to get the local physical density) 
		);

		converter.print();
		writeLogFile(converter, "phaseSeperation");

		//const T omega1 = 1.0;
		const T omega1 = converter.getOmega();
		const T G = -6.; //SC93
							  //const T G = -1/0.0649; //für t_sat = 0.0649 //PengRobinson


		// === 2rd Step: Prepare Geometry ===
		Vector<T, 2> extend(nx, ny);
		Vector<T, 2> origin(0, 0);
		IndicatorCuboid2D<T> cuboid(extend, origin);

		// Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
		const int noOfCuboids = singleton::mpi().getSize();
#else
		const int noOfCuboids = 1;
#endif
		CuboidGeometry2D<T> cuboidGeometry(cuboid, converter.getLatticeL(), noOfCuboids);
		//CuboidGeometry2D<T> cuboidGeometry(0, 0, 1, N, N, noOfCuboids); //hier int nx und int ny??

		// Periodic boundaries in x- and y-direction
		cuboidGeometry.setPeriodicity(true, true);

		// Instantiation of a loadBalancer
		HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

		// Instantiation of a superGeometry
		SuperGeometry2D<T> superGeometry(cuboidGeometry, loadBalancer, 2);

		prepareGeometry(converter, superGeometry, radius);

		// === 3rd Step: Prepare Lattice ===
		SuperLattice2D<T, DESCRIPTOR> sLattice(superGeometry);

		//VelocityShifting
		ForcedShanChenBGKdynamics<T, DESCRIPTOR> bulkDynamics1(
			omega1, instances::getExternalVelocityMomenta<T, DESCRIPTOR>());

		std::vector<T> rho0;
		rho0.push_back(1);
		rho0.push_back(1);
		//ShanChen93<T, T> interactionPotential;
		PengRobinson<T, T> interactionPotential(G);
		ShanChenForcedSingleComponentGenerator2D<T, DESCRIPTOR> coupling(G, rho0, interactionPotential);

		sLattice.addLatticeCoupling(superGeometry, 0, coupling, sLattice);
		sLattice.addLatticeCoupling(superGeometry, 1, coupling, sLattice);

		prepareLattice(converter, sLattice, bulkDynamics1, superGeometry);

		// === 4th Step: Main Loop ===

		//Konvergenzkriterium
		T interval = 20; //vergrößern
		T epsilon = 0.001;
		util::ValueTracer<T> converge(interval, epsilon);

		int iT = 0;

		clout << "starting simulation..." << endl;
		Timer<T> timer(maxIter, superGeometry.getStatistics().getNvoxel());
		timer.start();
		for (iT = 0; iT < maxIter; ++iT) {

			//Falls System konvergiert
			if (converge.hasConverged()) {
				clout << "Simulation converged." << endl;
				getResults(sLattice, converter, iT, superGeometry, timer, converge.hasConverged(), radius); //hier
				clout << "hier" << endl;
				break;
			}

			// === 5th Step: Definition of Initial and Boundary Conditions ===
			// in this application no boundary conditions have to be adjusted

			// === 6th Step: Collide and Stream Execution ===
			sLattice.collideAndStream();
			sLattice.communicate();
			sLattice.executeCoupling();

			// === 7th Step: Computation and Output of the Results ===
			getResults(sLattice, converter, iT, superGeometry, timer, converge.hasConverged(), radius);
			converge.takeValue(sLattice.getStatistics().getAverageEnergy(), true); //maxU
		}
		//hier dichte extrahieren und jeweiligen druck für oberflächenspannung berechnen
		timer.stop();
		timer.printSummary();
		iter++;
	}//while end
}



