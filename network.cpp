#include "network.hpp"

#include <fstream>
#include <cmath>
#include <boost/graph/graph_traits.hpp>

namespace dynamicNet
{	
	// graph definition
	graph Network::g;
	// number vertices/nodes
	std::size_t Network::N;
	// number edges
	std::size_t Network::edges;
	// degree
	double Network::degree;
	// network iterators
	graph::adjacency_iterator Network::neighbourIt, Network::neighbourEnd;
	graph::adjacency_iterator Network::nextNeighbourIt, Network::nextNeighbourEnd;

	Network::Network( std::string filepath )
	{
		
		std::ifstream inputFile(filepath.c_str());
		std::string skipline;
		std::getline(inputFile, skipline);
		
		while (true)
		{	
			int x, y;
			inputFile >> x >> y;
			if( inputFile.eof() ) break;
			boost::add_edge(x, y, g);
		}
		
		// number vertices/nodes
		N = boost::num_vertices(g);
		// number edges
		edges = boost::num_edges(g);
		// degree
		degree = 1.0*edges/N;
	}

} // end namespace dynamicNet
