#ifndef KMC_HPP
#define KMC_HPP

#include <random>
#include <vector>

namespace dynamicNet
{

extern std::mt19937 rng;
	
// kinetic Monte Carlo
class kineticMC
{
public:
	
	// kmc variables
	static double kmcTime_,Z_,IntZR_,ExtZR_,IntFR_,ExtFR_;
	static unsigned mcount;
		
	// initial value functions
	static void set_internal_failure_rate( double );
	static void set_external_failure_rate( double );
	static void set_internal_recovery_rate( double );
	static void set_external_recovery_rate( double );
	static void set_neighborhood_failure_value( unsigned );

	// hash lists to show the current node state
	static std::vector<unsigned> failureHashList_;
	static std::vector<int> externalFailureHashList_;
	static std::vector<int> internalFailureHashList_;

	// rates
	static double R_,Q1_,Q2_,Q3_,Q4_;
	
	// place ring
	static void place_ring( unsigned );
	
	// simulation step
	static void sim_step();

	// default constructor
	kineticMC();
	
private:
	
	// shuffle rng
	static int shuffle_mtrng (int n);
	
	// update indices
	static void update_indices( std::vector<int>&, std::vector<unsigned>&, unsigned );

	// simulation functions
	static void internalFailure( unsigned );
	static void externalFailure( unsigned );
	static void internalRecovery( unsigned );
	static void externalRecovery( unsigned );

	// input values
	static double internal_failure_rate_;
	static double external_failure_rate_;
	static double internal_recovery_rate_;
	static double external_recovery_rate_;
	static unsigned neighborhood_failure_value_;
	
	// simulation rng
	static std::uniform_real_distribution<double> kmcRand;	
	
	// vector containers
	static std::vector<unsigned> internalFailureRateList_;
	static std::vector<unsigned> externalFailureRateList_;
	static std::vector<unsigned> internalActiveRateList_;
	static std::vector<unsigned> externalActiveRateList_;


};
} // end namespace dynamicNet

#endif // !defined KMC_HPP
