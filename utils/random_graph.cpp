#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/random.hpp>
#include <iostream>
#include <fstream>
#include <boost/graph/graphviz.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

int main()
{
    using namespace boost;
    using namespace boost::random;
    using namespace std;
    using Graph = adjacency_list<vecS, vecS, directedS>;

     // Define the parameters for the random graph
    int num_vertices = 10;
    int num_edges = 15;

    // Seed the random number generator
    boost::random::mt19937 rng; 
    rng.seed(time(0));

    Graph random_graph;
    boost::generate_random_graph(random_graph, num_vertices, num_edges, rng);

    std::ofstream dotfile("TestNetwork.dot");
    write_graphviz(dotfile, random_graph);
    std::cout << "Graph visualization exported to TestNetwork.dot" << std::endl;
    
    return 0;
}