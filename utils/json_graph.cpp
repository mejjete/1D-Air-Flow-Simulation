#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/graphviz.hpp>
#include <tuple>
#include <algorithm>
#include "network.hpp"

int main(int argc, char **argv)
{
    using namespace boost;
    using namespace std;

    // User must specify JSON configuration file for entire network
    if(argc <= 1)
    {
        std::cerr << "No input file exist\n";
        exit(EXIT_FAILURE);
    }

    boost::property_tree::ptree json_network;
    boost::property_tree::read_json(argv[1], json_network);

    const std::string network_name = json_network.get<std::string>("name");
    const int dx = json_network.get<int>("dx");
    const double h = json_network.get<double>("h");
    const int t_step = json_network.get<int>("t_step");
    const int total_vert = json_network.get<int>("vertices");
    const int total_edg = json_network.get<int>("edges");

    auto default_out = cout.flags();
    cout << "--------------------------------------------------" << endl;
    cout << "Simulation parameters" << endl;
    cout << "--------------------------------------------------" << endl;
    cout.flags(std::ios::left);
    cout << setw(10) << "dx:" << setw(10) << dx << endl;
    cout << setw(10) << "dt:" << setw(10) << h << endl;
    cout << setw(10) << "time:" << setw(10) << t_step << endl;
    cout << "--------------------------------------------------" << endl;
    cout.flags(default_out);

    /**
     * Create graph with specified number of vertices, then connect them 
     * according to a JSON scheme.
     */
    using Graph = adjacency_list<vecS, vecS, directedS, VertexProperty, EdgeProperty>;
    using EdgePropertyMap = typename property_map<Graph, edge_bundle_t>::type;
    using VertexPropertyMap = typename property_map<Graph, vertex_bundle_t>::type;
    Graph network;

    std::vector<std::pair<int, double>> init_pressure;
    boost::property_tree::ptree arrayPressure = json_network.get_child("vens");

    // Fetch number of vens their pressure, and associated vertex ID
    for(auto &child : arrayPressure)
    {
        int vertex_ID = child.second.get<int>("vertex_id");
        double pressure = child.second.get<double>("pressure");
        init_pressure.push_back(std::make_pair(vertex_ID, pressure));
    }

    for(int i = 0; i < total_vert; i++)
    {
        /**
         *  Initial pressure for a single edge is stored inside it's targer vertex. 
         *  Following code sets up initial pressure for vertices where we have vens. 
         *  For other edges pressure set to 0 because it will be evaluated during simulation.
        */
        auto is_init_P = [i](std::pair<int, double> item) 
        {
            if(i == item.first)
                return true;
            return false;
        };

        auto item = find_if(init_pressure.begin(), init_pressure.end(), is_init_P);
        if(item != init_pressure.end())
            add_vertex(VertexProperty(i, h, item->second), network);
        else 
            add_vertex(VertexProperty(i, h, 0.0), network);
    }
    
    // Initialize each edge and connect it to the vertex
    int total_edges = 0;
    
    boost::property_tree::ptree arrayEdges = json_network.get_child("edge");
    for(auto &child : arrayEdges)
    {
        int id = child.second.get<int>("id");
        int length = (int)(ceil(child.second.get<double>("length")));
        double alpha = child.second.get<double>("alpha");
        double beta = child.second.get<double>("beta");
        double gamma = child.second.get<double>("gamma");
        int s_step = (int)(ceil((double)length / (double)dx));
        int head = child.second.get<int>("head");
        int tail = child.second.get<int>("tail");

        EdgeProperty insert_edge(id, length, s_step, t_step, alpha, beta, gamma, h);
        add_edge(head, tail, insert_edge, network);
        total_edges++;
    }

    using EdgeIter = typename boost::graph_traits<Graph>::edge_iterator;
    using VertexIter = typename boost::graph_traits<Graph>::vertex_iterator;

    EdgePropertyMap edge_map = get(edge_bundle, network);
    VertexPropertyMap vertex_map = get(vertex_bundle, network);

    /**
     *  Each vertex has some number of incoming edges. For a given vertex, set the gamma that 
     *  is biggest among incoming edges. 
    */
    for(std::pair<EdgeIter, EdgeIter> iter = edges(network); iter.first != iter.second; ++iter.first)
    {
        auto tar_vert = target(*iter.first, network);
        EdgeProperty &edge = network[*iter.first];
        VertexProperty &vertex = network[tar_vert];

        if(edge.getGamma() > vertex.getGamma())
            vertex.setGamma(edge.getGamma());
    }

    if(total_edges != total_edg)
    {
        std::cerr << "Total number of edges does not match\n";
        exit(EXIT_FAILURE);
    }

    std::ofstream dotfile("TestNetwork.dot");
    write_graphviz(dotfile, network);
    std::cout << "Graph visualization exported to TestNetwork.dot" << std::endl;

    return 0;
}