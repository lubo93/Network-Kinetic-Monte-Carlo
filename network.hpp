#ifndef NETWORK_HPP
#define NETWORK_HPP

#include <iostream>
#include <boost/graph/adjacency_list.hpp>

namespace dynamicNet
{
	
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> graph; 

	// network
	class Network
	{
	public:
		
		// graph definition
		static graph g;
		// number vertices/nodes
		static std::size_t N;
		// number edges
		static std::size_t edges;
		// network degree
		static double degree;
		// network iterators
		static graph::adjacency_iterator neighbourIt, neighbourEnd;
		static graph::adjacency_iterator nextNeighbourIt, nextNeighbourEnd;

		// Default constructor
		Network( std::string );
		
	};
} // end namespace dynamicNet

#endif // !defined NETWORK_HPP
