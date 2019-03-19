// dynamicNet Simulation
// Lucas Boettcher
// http://www.lucas-boettcher.info

#include "network.hpp"
#include "kmc.hpp"

#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <cstdlib>
#include <random>
#include <ctime>

#ifdef _MPI
	#include "mpi.h"
#endif

using namespace dynamicNet;

int main(int argc, char** argv)
{	
	
#ifdef _MPI
	int id;
	int num;
	//
	//  initialize MPI
	//
	MPI_Init ( &argc, &argv );
	//
	// processes number
	//
	MPI_Comm_size(MPI_COMM_WORLD, &num);
	//
	//  individual process ID
	//
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	//
#endif
	
	if (argc!=9) {
        std::cout << argc << "Specify input arguments after " << argv[0] << " : " << "[time T] [measurement steps nmeas]"
		<< "[internal failure rate p] [external failure rate r] [neighborhood failure value m] "
		<< "[internal recovery rate QInt]" "[external recovery rate QExt]" << std::endl;
        return -1;
    }
    std::size_t sim_time=std::stoul(argv[1]),nmeas=std::stoul(argv[2]);
	double p=std::stod(argv[3]),r=std::stod(argv[4]),m=std::stod(argv[5]),qint=std::stod(argv[6]),qext=std::stod(argv[7]);
	std::string filepath=argv[8];
	
	kineticMC::set_internal_failure_rate		( p );
	kineticMC::set_external_failure_rate		( r );
	kineticMC::set_neighborhood_failure_value	( m );
	kineticMC::set_internal_recovery_rate		( qint );
	kineticMC::set_external_recovery_rate		( qext );


	// define network
	Network dynamicNet( filepath );
	// define kinetic Monte Carlo
	// for given network
	kineticMC kmcSIS;

#ifdef _MPI	
	char out_filename[160];
	sprintf( out_filename, "dynamicNetdata_N%d_p%05.3f_r%05.3f_%04d.dat", unsigned(dynamicNet.N), p, r, id );
	std::ofstream outfile(out_filename);
#else
	char out_filename[160];
	sprintf( out_filename, "dynamicNetdata_N%d_p%05.3f_r%05.3f.dat", unsigned(dynamicNet.N), p, r );
	std::ofstream outfile(out_filename);
#endif
	
// random seed for production code
#ifdef _NDEBUG
	auto const seed = std::random_device()();
    dynamicNet::rng.seed(seed);
	outfile << "# Seed: " << seed << std::endl;
#endif
	
	outfile << "# N: " << dynamicNet.N << std::endl;
	outfile << "# p: " << p << std::endl;
	outfile << "# r: " << r << std::endl;
	outfile << "# m: " << m << std::endl;
	outfile << "# qint: " << qint << std::endl;
	outfile << "# qext: " << qext << std::endl;
	outfile << "# k: " << dynamicNet.degree << std::endl;
	outfile << "#" << std::endl;

	kmcSIS.place_ring(0);
	
	unsigned cnt = 0;
	clock_t begin = clock();
	while ( kmcSIS.kmcTime_ <= sim_time)
	{	
		if (cnt%nmeas==0)
			outfile << kmcSIS.kmcTime_ << " " << kmcSIS.Z_/dynamicNet.N << " " << kmcSIS.ExtZR_/dynamicNet.N << " " << kmcSIS.IntZR_/dynamicNet.N  << " " << kmcSIS.ExtFR_/dynamicNet.N << " " << kmcSIS.IntFR_/dynamicNet.N << std::endl;

		// KMC simulation step
		kmcSIS.sim_step();
		
		++cnt;
	}
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	outfile << "#"  << std::endl;
	outfile << "# Time needed: " << elapsed_secs << std::endl;
	std::cout << "Time needed: " << elapsed_secs << std::endl;

#ifdef _MPI
	//
	//  terminate MPI
	//
	MPI_Finalize ( );
	//
	//  terminated MPI
	//
#endif
	
	return 0;
}
