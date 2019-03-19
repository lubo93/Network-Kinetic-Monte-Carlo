#include <cassert>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "kmc.hpp"
#include "network.hpp"

namespace dynamicNet
{
	// rng initialization
	std::mt19937 rng(0);
	
	// input values
	double kineticMC::internal_failure_rate_ = 0;
	double kineticMC::external_failure_rate_ = 0;
	double kineticMC::internal_recovery_rate_ = 0;
	double kineticMC::external_recovery_rate_ = 0;
	unsigned kineticMC::neighborhood_failure_value_ = 0;

	// model rates (initialized to 1)
	// because of main loop
	double kineticMC::R_ = 1;
	double kineticMC::Q1_ = 1;
	double kineticMC::Q2_ = 1;
	double kineticMC::Q3_ = 1;
	double kineticMC::Q4_ = 1;
	
	// global variables
	double kineticMC::kmcTime_ = 0;
	double kineticMC::Z_ = 0;
	double kineticMC::IntZR_ = 0;
	double kineticMC::ExtZR_ = 0;
	double kineticMC::IntFR_ = 0;
	double kineticMC::ExtFR_ = 0;
	unsigned kineticMC::mcount = 0;

	// vector containers
	std::vector<unsigned> kineticMC::internalFailureRateList_;
	std::vector<unsigned> kineticMC::externalFailureRateList_;
	std::vector<unsigned> kineticMC::internalActiveRateList_;
	std::vector<unsigned> kineticMC::externalActiveRateList_;
	
	// hash lists
	// failure hash list: 0: active and 1: failure
	// external failure hash list: -1: no external failure
	// and >-1: index in externalActiveRateList_
	std::vector<unsigned> kineticMC::failureHashList_;
	std::vector<int> kineticMC::externalFailureHashList_;
	std::vector<int> kineticMC::internalFailureHashList_;

	// simulation rng
	std::uniform_real_distribution<> kineticMC::kmcRand(0,1);
	
	// RNG for vector shuffling
	int kineticMC::shuffle_mtrng (int n) 
	{	
		std::uniform_int_distribution<unsigned> initial_distribution(0,n-1);
		return initial_distribution(rng);
	}
	
	// set initial values
	void kineticMC::set_internal_failure_rate( double p )
	{
		internal_failure_rate_ = p;
	}

	void kineticMC::set_external_failure_rate( double r )
	{
		external_failure_rate_ = r;
	}

	void kineticMC::set_neighborhood_failure_value( unsigned m )
	{
		neighborhood_failure_value_ = m;
	}
	
	void kineticMC::set_internal_recovery_rate( double qint )
	{
		internal_recovery_rate_ = qint;
	}

	void kineticMC::set_external_recovery_rate( double qext )
	{
		external_recovery_rate_ = qext;
	}
	
	// update indices after element has been erased
	void kineticMC::update_indices( std::vector<int> &v, std::vector<unsigned> &w, unsigned erased )
	{
		unsigned node_back = w.back();
		unsigned node_erased = w[erased];
		
		w[erased] = node_back;
		w.pop_back();
		v[node_back] = erased;
		v[node_erased] = -1;
	}

	// internal failure
	void kineticMC::internalFailure( unsigned node )
	{
		// update failure hash lists
		failureHashList_[node] = 1;
		internalActiveRateList_.pop_back();
		internalFailureHashList_[node] = -1;
		internalFailureRateList_.push_back(node);

		// update list numbers
		--Z_;
		--IntZR_;
		++IntFR_;

		// for this node check if it was exposed
		// to external failure (number active neighbours <= m)
		if (externalFailureHashList_[node] > -1)
		{
			// since node erased update
			// hash list nodes with index > node index
			update_indices( externalFailureHashList_, externalActiveRateList_, externalFailureHashList_[node] );
			--ExtZR_;
		}

		// for all neighbours of this node check
		// if they are now exposed to external failure
		// (number active neighbours <= m)
		boost::tie(Network::neighbourIt, Network::neighbourEnd) = boost::adjacent_vertices( node, Network::g );
		for (; Network::neighbourIt != Network::neighbourEnd; ++Network::neighbourIt)
		{
			mcount = 0;
			if (failureHashList_[*Network::neighbourIt] == 0 && externalFailureHashList_[*Network::neighbourIt] == -1)
			{
				boost::tie(Network::nextNeighbourIt, Network::nextNeighbourEnd) = boost::adjacent_vertices( *Network::neighbourIt, Network::g );
				for (; Network::nextNeighbourIt != Network::nextNeighbourEnd; ++Network::nextNeighbourIt)
				{
					if (failureHashList_[*Network::nextNeighbourIt] == 0)	++mcount;
				}
				if (mcount <= neighborhood_failure_value_)
				{
					externalActiveRateList_.push_back(*Network::neighbourIt);
					externalFailureHashList_[*Network::neighbourIt] = externalActiveRateList_.size()-1;
					++ExtZR_;
				}
			}
		}

	}

	// external failure
	void kineticMC::externalFailure( unsigned node )
	{
		// update failure hash lists
		failureHashList_[node] = 1;
		externalActiveRateList_.pop_back();
		externalFailureHashList_[node] = -1;
		externalFailureRateList_.push_back(node);
		// since node erased update
		// hash list nodes with index > node index
		update_indices( internalFailureHashList_, internalActiveRateList_, internalFailureHashList_[node]);

		// update list numbers
		--Z_;
		--ExtZR_;
		--IntZR_;
		++ExtFR_;

		// for all neighbours of this node check
		// if they are now exposed to external failure
		// (number active neighbours <= m)
		boost::tie(Network::neighbourIt, Network::neighbourEnd) = boost::adjacent_vertices( node, Network::g );
		for (; Network::neighbourIt != Network::neighbourEnd; ++Network::neighbourIt)
		{
			mcount = 0;
			if (failureHashList_[*Network::neighbourIt] == 0 && externalFailureHashList_[*Network::neighbourIt] == -1)
			{
				boost::tie(Network::nextNeighbourIt, Network::nextNeighbourEnd) = boost::adjacent_vertices( *Network::neighbourIt, Network::g );
				for (; Network::nextNeighbourIt != Network::nextNeighbourEnd; ++Network::nextNeighbourIt)
				{
					if (failureHashList_[*Network::nextNeighbourIt] == 0)	++mcount;
				}
				if (mcount <= neighborhood_failure_value_)
				{
					externalActiveRateList_.push_back(*Network::neighbourIt);
					externalFailureHashList_[*Network::neighbourIt] = externalActiveRateList_.size()-1;
					++ExtZR_;
				}
			}
		}
	}

	// internal recovery
	void kineticMC::internalRecovery( unsigned node )
	{

		// update failure hash lists
		failureHashList_[node] = 0;
		internalFailureRateList_.pop_back();
		internalActiveRateList_.push_back(node);
		internalFailureHashList_[node] = internalActiveRateList_.size()-1;

		// update list numbers
		++Z_;
		++IntZR_;
		--IntFR_;

		// for this node check if it is exposed
		// to external failure (number active neighbours <= m)
		mcount = 0;
		boost::tie(Network::neighbourIt, Network::neighbourEnd) = boost::adjacent_vertices( node, Network::g );
		for (; Network::neighbourIt != Network::neighbourEnd; ++Network::neighbourIt)
		{
			if (failureHashList_[*Network::neighbourIt] == 0)	++mcount;
		}
		if (mcount <= neighborhood_failure_value_)
		{
			externalActiveRateList_.push_back(node);
			externalFailureHashList_[node] = externalActiveRateList_.size()-1;
			++ExtZR_;
		}

		// for all neighbours of this node check
		// if they are now not exposed to external failure anymore
		// (number active neighbours > m)
		boost::tie(Network::neighbourIt, Network::neighbourEnd) = boost::adjacent_vertices( node, Network::g );
		for (; Network::neighbourIt != Network::neighbourEnd; ++Network::neighbourIt)
		{
			mcount = 0;
			if (failureHashList_[*Network::neighbourIt] == 0 && externalFailureHashList_[*Network::neighbourIt] > -1)
			{
				boost::tie(Network::nextNeighbourIt, Network::nextNeighbourEnd) = boost::adjacent_vertices( *Network::neighbourIt, Network::g );
				for (; Network::nextNeighbourIt != Network::nextNeighbourEnd; ++Network::nextNeighbourIt)
				{
					if (failureHashList_[*Network::nextNeighbourIt] == 0)	++mcount;
				}
				if (mcount > neighborhood_failure_value_)
				{
					// since node erased update
					// hash list nodes with index > node index
					update_indices( externalFailureHashList_, externalActiveRateList_, externalFailureHashList_[*Network::neighbourIt] );
					--ExtZR_;
				}
			}
		}
	}

	// external recovery
	void kineticMC::externalRecovery( unsigned node )
	{
		// update failure hash lists
		failureHashList_[node] = 0;
		externalFailureRateList_.pop_back();
		internalActiveRateList_.push_back(node);
		internalFailureHashList_[node] = internalActiveRateList_.size()-1;

		// update list numbers
		++Z_;
		++IntZR_;
		--ExtFR_;

		// for this node check if it is exposed
		// to external failure (number active neighbours <= m)
		mcount = 0;
		boost::tie(Network::neighbourIt, Network::neighbourEnd) = boost::adjacent_vertices( node, Network::g );
		for (; Network::neighbourIt != Network::neighbourEnd; ++Network::neighbourIt)
		{
			if (failureHashList_[*Network::neighbourIt] == 0)	++mcount;
		}
		if (mcount <= neighborhood_failure_value_)
		{
			externalActiveRateList_.push_back(node);
			externalFailureHashList_[node] = externalActiveRateList_.size()-1;
			++ExtZR_;
		}

		// for all neighbours of this node check
		// if they are now not exposed to external failure anymore
		// (number active neighbours > m)
		boost::tie(Network::neighbourIt, Network::neighbourEnd) = boost::adjacent_vertices( node, Network::g );
		for (; Network::neighbourIt != Network::neighbourEnd; ++Network::neighbourIt)
		{
			mcount = 0;
			if (failureHashList_[*Network::neighbourIt] == 0 && externalFailureHashList_[*Network::neighbourIt] > -1)
			{
				boost::tie(Network::nextNeighbourIt, Network::nextNeighbourEnd) = boost::adjacent_vertices( *Network::neighbourIt, Network::g );
				for (; Network::nextNeighbourIt != Network::nextNeighbourEnd; ++Network::nextNeighbourIt)
				{
					if (failureHashList_[*Network::nextNeighbourIt] == 0)	++mcount;
				}
				if (mcount > neighborhood_failure_value_)
				{
					// since node erased update
					// hash list nodes with index > node index
					update_indices( externalFailureHashList_, externalActiveRateList_, externalFailureHashList_[*Network::neighbourIt] );
					--ExtZR_;
				}
			}
		}
	}

	// simulation step
	void kineticMC::sim_step()
	{
		// determine the rates
		Q1_ = IntZR_*internal_failure_rate_;
		Q2_ = ExtZR_*external_failure_rate_;
		Q3_ = IntFR_*internal_recovery_rate_;
		Q4_ = ExtFR_*external_recovery_rate_;
		R_ = Q1_+Q2_+Q3_+Q4_;
		//std::cout << kmcTime_ << " " << Z_ << " " << internalActiveRateList_.size() << " " << externalActiveRateList_.size() << " " << internalFailureRateList_.size() << " " << externalFailureRateList_.size() << std::endl;
		//std::cout << IntZR_ << " " << ExtZR_ << " " << IntFR_ << " " << ExtFR_ << " " << std::endl;
		
		// failure or recovery?
		if (kmcRand(rng) < (Q1_+Q2_)/R_)
		{
			// internal or external failure?
			if (kmcRand(rng) < Q1_/(Q1_+Q2_))
			{
				unsigned random_index = shuffle_mtrng(internalActiveRateList_.size());
				unsigned node_random = internalActiveRateList_[random_index];
				unsigned node_back = internalActiveRateList_.back();
				
				internalActiveRateList_[random_index] = node_back;
				internalActiveRateList_.pop_back();
				internalActiveRateList_.push_back(node_random);
				
				internalFailureHashList_[node_back] = random_index;
				internalFailureHashList_[node_random] = internalActiveRateList_.size()-1;;

				internalFailure( internalActiveRateList_.back() );
			}
			else
			{	
				unsigned random_index = shuffle_mtrng(externalActiveRateList_.size());
				unsigned node_random = externalActiveRateList_[random_index];
				unsigned node_back = externalActiveRateList_.back();
				
				externalActiveRateList_[random_index] = node_back;
				externalActiveRateList_.pop_back();
				externalActiveRateList_.push_back(node_random);
				
				externalFailureHashList_[node_back] = random_index;
				externalFailureHashList_[node_random] = externalFailureHashList_.size()-1;
			
				externalFailure( externalActiveRateList_.back() );
			}
		}
		else
		{
			// internal or external recovery?
			if (kmcRand(rng) < Q3_/(Q3_+Q4_))
			{
				unsigned random_index = shuffle_mtrng(internalFailureRateList_.size());
				unsigned node_random = internalFailureRateList_[random_index];
				unsigned node_back = internalFailureRateList_.back();
				
				internalFailureRateList_[random_index] = node_back;
				internalFailureRateList_.pop_back();
				internalFailureRateList_.push_back(node_random);
			
				internalRecovery( internalFailureRateList_.back() );
			}
			else
			{
				unsigned random_index = shuffle_mtrng(externalFailureRateList_.size());
				unsigned node_random = externalFailureRateList_[random_index];
				unsigned node_back = externalFailureRateList_.back();
				
				externalFailureRateList_[random_index] = node_back;
				externalFailureRateList_.pop_back();
				externalFailureRateList_.push_back(node_random);
				
				externalRecovery( externalFailureRateList_.back() );
			}
		}

		kmcTime_ += -1.0/R_*log(1-kmcRand(rng));

	}
	
	// simulation step
	void kineticMC::place_ring( unsigned radius )
	{
		// spreading ring construction
		std::vector<unsigned> spreadingRing;
		unsigned center = int(Network::N/2+sqrt(Network::N)/2);
		std::cout << center << std::endl;
		// define spreading ring
		//for (std::size_t i = 0; i < (2*radius+1); ++i)
		//{
		//    for (std::size_t j = 0; j < (2*radius+1); ++j)
		//    {
		//	if (pow(abs(radius-i),2)+pow(abs(radius-j),2)-pow(radius,2) <= 0)
		//	{
		//	  spreadingRing.push_back(center-(radius-i)*int(sqrt(Network::N))-radius+j);
		//	}
		//    }
		//}
		
		std::cout << 3.14*pow(radius,2)/Network::N << std::endl;
		for (std::size_t i = 0; i < radius; ++i)
		{
			spreadingRing.push_back(i);
		}
		
		// update nodes inside spreading ring
		for (std::size_t i = 0; i < spreadingRing.size(); ++i)
		{
			// failure inside spreading ring is externally caused
			// check if node is already infected
			if (failureHashList_[spreadingRing[i]] == 1)
			{
			    // change internal to external failure
			    for (std::size_t j = 0; j < internalFailureRateList_.size(); ++j)
			    {
				if (spreadingRing[i] == internalFailureRateList_[j])
				{
				     unsigned node_back = internalFailureRateList_.back();
				     internalFailureRateList_[j] = node_back;
				     internalFailureRateList_.pop_back();
				     externalFailureRateList_.push_back(spreadingRing[i]);
				     --IntFR_;
				     ++ExtFR_;
				}
			    }
			}
			else
			{
			    // set failure (assume externally caused)
			    failureHashList_[spreadingRing[i]] = 1;
			    // modify internal lists
			    unsigned node_back = internalActiveRateList_.back();
			    internalActiveRateList_.pop_back();
			    internalActiveRateList_[internalFailureHashList_[spreadingRing[i]]] = node_back;
			    internalFailureHashList_[node_back] = internalFailureHashList_[spreadingRing[i]];
			    internalFailureHashList_[spreadingRing[i]] = -1;
			    
			    // let the node externally fail
			    externalFailureRateList_.push_back(spreadingRing[i]);
			    
			    // update list numbers
			    --Z_;
			    --IntZR_;
			    ++ExtFR_;
			    
			    // for this node check if it was exposed
			    // to external failure (number active neighbours <= m)
			    if (externalFailureHashList_[spreadingRing[i]] > -1)
			    {
				    // since node erased update
				    // hash list nodes with index > node index
				    update_indices( externalFailureHashList_, externalActiveRateList_, externalFailureHashList_[spreadingRing[i]] );
				    --ExtZR_;
			    }
			    
			}
		}
		
		// update nodes inside spreading ring
		for (std::size_t i = 0; i < spreadingRing.size(); ++i)
		{
			// for all neighbours of this node check
			// if they are now exposed to external failure
			// (number active neighbours <= m)
			boost::tie(Network::neighbourIt, Network::neighbourEnd) = boost::adjacent_vertices( spreadingRing[i], Network::g );
			for (; Network::neighbourIt != Network::neighbourEnd; ++Network::neighbourIt)
			{
				mcount = 0;
				if (failureHashList_[*Network::neighbourIt] == 0 && externalFailureHashList_[*Network::neighbourIt] == -1)
				{
					boost::tie(Network::nextNeighbourIt, Network::nextNeighbourEnd) = boost::adjacent_vertices( *Network::neighbourIt, Network::g );
					for (; Network::nextNeighbourIt != Network::nextNeighbourEnd; ++Network::nextNeighbourIt)
					{
						if (failureHashList_[*Network::nextNeighbourIt] == 0)	++mcount;
					}
					if (mcount <= neighborhood_failure_value_)
					{
						externalActiveRateList_.push_back(*Network::neighbourIt);
						externalFailureHashList_[*Network::neighbourIt] = externalActiveRateList_.size()-1;
						++ExtZR_;
					}
				}
			}
		}
	}
	
	// default constructor
	kineticMC::kineticMC()
	{	
	  	kmcTime_ = 0;
		Z_ = 0;
		IntZR_ = 0;
		ExtZR_ = 0;
		IntFR_ = 0;
		ExtFR_ = 0;
		mcount = 0;
	
		internalFailureRateList_.clear();
		externalFailureRateList_.clear();
		internalActiveRateList_.clear();
		externalActiveRateList_.clear();
		
		failureHashList_.clear();
		externalFailureHashList_.clear();
		internalFailureHashList_.clear();
	
		// initialize containers
		for (std::size_t i = 0; i < Network::N; ++i)
		{
			failureHashList_.push_back(0);
			externalFailureHashList_.push_back(-1);
			internalActiveRateList_.push_back(i);
			internalFailureHashList_.push_back(i);
			++IntZR_;
			++Z_;
		}
		
		// initial values
		std::cout << std::endl;
		std::cout << "Initial number of active nodes (Z): " << Z_ << std::endl;
		std::cout << "Initial internal active rate (IntZR): " << IntZR_ << std::endl;
		std::cout << "Initial external active rate (ExtZR): " << ExtZR_ << std::endl;
		std::cout << std::endl;


	}
} // end namespace dynamicNet
