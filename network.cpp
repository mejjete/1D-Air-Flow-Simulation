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
    const double rho = json_network.get<double>("rho");
    const int a = json_network.get<int>("a");

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
    using Graph = adjacency_list<vecS, vecS, directedS, VertexProperty*, EdgeProperty*>;
    using EdgePropertyMap = typename property_map<Graph, edge_bundle_t>::type;
    using VertexPropertyMap = typename property_map<Graph, vertex_bundle_t>::type;
    Graph network;

    std::vector<std::pair<int, double>> init_pressure;
    boost::property_tree::ptree arrayPressure = json_network.get_child("vens");

    // Fetch number of vens, their pressure and associated vertex ID
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
        {
            VertexProperty *insert = new VertexProperty(i, h, item->second); 
            add_vertex(insert, network);
        }
        else 
        {
            VertexProperty *insert = new VertexProperty(i, h, 0.0); 
            add_vertex(insert, network);
        }
    }
    
    // Initialize each edge and connect it to the vertex
    int total_edges = 0;
    
    boost::property_tree::ptree arrayEdges = json_network.get_child("edge");
    for(auto &child : arrayEdges)
    {
        int id = child.second.get<int>("id");
        int length = (int)(ceil(child.second.get<double>("length")));
        double F = child.second.get<double>("F");
        double r = child.second.get<double>("r");
        double alpha = F / (rho * dx);
        double beta = (F * r) / rho; 
        double gamma = (rho * a * a) / (F * dx);
        int s_step = (int)(ceil((double)length / (double)dx));
        int head = child.second.get<int>("head");
        int tail = child.second.get<int>("tail");

        EdgeProperty *insert_edge = new EdgeProperty(id, s_step, t_step, alpha, beta, gamma, h);
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
        EdgeProperty *edge = network[*iter.first];
        VertexProperty *vertex = network[tar_vert];

        if(edge->getGamma() > vertex->getGamma())
            vertex->setGamma(edge->getGamma());
    }

    if(total_edges != total_edg)
    {
        std::cerr << "Total number of edges does not match\n";
        exit(EXIT_FAILURE);
    }

    // Export graph to a Graphviz DOT file
    #ifdef DEBUG
    std::map<Graph::edge_descriptor, int> edge_numbers;

    // Assign numbers to each edge
    int edge_num = 0;
    for (auto e : boost::make_iterator_range(boost::edges(network))) {
        edge_numbers[e] = edge_num++;
    }

    std::ofstream dotfile("TestNetwork.dot");
    write_graphviz(dotfile, network);
    std::cout << "Graph visualization exported to TestNetwork.dot" << std::endl;
    
    std::ofstream edgefile("TestNetwork_Vertices.txt");
    edgefile << "Vertex set: " << std::endl;

    for(auto iter = vertices(network); iter.first != iter.second; ++iter.first)
    {
        VertexProperty *vertex = vertex_map[*iter.first];
        edgefile << "Gamma: " << vertex->getGamma() << endl;
    }
    std::cout << "Vertex information exported to TestNetwork_Vertices.txt" << std::endl;
    #endif 


    /**
     *  Calculation of air flow and pressure for entire network
     * 
     * 
     * 
     *  Phase 1: Flow and pressure calculation inside edges. 
     *  During this phase we calculate two parameters for each edge such as an independent entity.                                                             
     * 
     *  Phase 2: Pressure adaptation. 
     *  After initial calculations is done, we have to adapt pressure in connecting points 
     *  (vertices) according to pressure in neighboring points (vertices).
    */


    #ifdef DEBUG
    std::vector<EdgeDebug> edge_debug;
    std::vector<VertexDebug> vertex_debug;

    auto edge_finder = [&edge_debug](int id) -> typename std::add_lvalue_reference<EdgeDebug>::type
    {
        for(auto &iter : edge_debug)
        {
            if(iter.getID() == id)
                return iter;
        }

        throw std::runtime_error("Edge debug array is corrupted\n");
    };

    auto vertex_finder = [&vertex_debug](int id) -> typename std::add_lvalue_reference<VertexDebug>::type
    {
        for(auto &iter : vertex_debug)
        {
            if(iter.getID() == id)
                return iter;
        }

        throw std::runtime_error("Vertex debug array is corrupted\n");
    };
    #endif

    for (int i = 0; i < t_step; i++) 
    {
        for(std::pair<EdgeIter, EdgeIter> et = edges(network); et.first != et.second; ++et.first)
        {
            auto tar_vrt_desc = target(*et.first, network);
            auto src_vrt_desc = source(*et.first, network);
            
            double target_pressure = network[tar_vrt_desc]->getP(i);
            double source_pressure = network[src_vrt_desc]->getP(i);

            /**
             *  Pressure initialization occurs at every timestep for each edge.
             *  Edges take pressure from source vertex at the start of simulation
             *  and commits last pressure to a target vertex when all calculations for edge 
             *  is done. This is the way they communicate via vertex pressure.
            */
            EdgeProperty *edge = edge_map[*et.first];
            edge->setSourceP(i, source_pressure);
            edge->setTargetP(i, target_pressure);

            int s_step = edge->getSteps();

            /*************************************************************/

            // Loop over spatial steps for flow
            for(int k = 1; k < s_step; k++)
                edge->calculateQ(i, k);

            // Loop over spatial steps for pressure
            for(int k = 1; k < s_step - 1; k++)
                edge->calculateP(i, k);

            /*************************************************************/

            VertexProperty *target_vertex = vertex_map(tar_vrt_desc);
            VertexProperty *source_vertex = vertex_map(src_vrt_desc);

            auto last_Q = edge->getLastQ(i);

            /**
             *  Pressure at first spatial step goes to a source vertex with a minus sign because it is 
             *  outcoming pressure relative to a source vertex.
             *  Pressure at last spatial step goes to a target vertex as a positive value because it is
             *  incoming pressure relative to a target vertex.
             * 
             *  See Kirchhoff's 1st law.
            */
            double first_flow = last_Q[0];
            double last_flow = last_Q[2];

            source_vertex->addQ(-first_flow);
            target_vertex->addQ(last_flow);

            #ifdef DEBUG
            // Initialize edge debug handlers
            if(i == 0)
                edge_debug.push_back(EdgeDebug("E" + std::to_string(edge->getID()), edge));

            // Write down time points for each edge separately
            if((i % (t_step / 3840)) == 0)
            {
                auto &current_debug = edge_finder(edge->getID());
                current_debug.serialize(i * h);
            }
            #endif
        }

        // Pressure adaptation in vertices
        for(auto vt = vertices(network); vt.first != vt.second; ++vt.first)
        {
            VertexProperty *vert = vertex_map[*vt.first];

            auto is_init_P = [vert](std::pair<int, double> item) 
            {
                if(vert->getID() == item.first)
                    return true;
                return false;
            };
            
            auto iter = find_if(init_pressure.begin(), init_pressure.end(), is_init_P);
            
            // Do not perform adaptation for vertices that have vans on it
            if(iter != init_pressure.end())
                continue; 
            vert->adapt(i);

            #ifdef DEBUG
            // Initialize vertex debug handlers
            if(i == 0)
                vertex_debug.push_back(VertexDebug("V" + std::to_string(vert->getID()), vert));

            // Write down time points for each vertex separately
            if((i % (t_step / 3840)) == 0)
            {
                auto &curr_vert_debug = vertex_finder(vert->getID());
                curr_vert_debug.serialize(i * h);
            }
            #endif
        }
    }

    for(std::pair<EdgeIter, EdgeIter> et = edges(network); et.first != et.second; ++et.first)
    {
        EdgeProperty *edge = edge_map[*et.first];
        auto tar_vrt_desc = target(*et.first, network);
        VertexProperty *vert = vertex_map[tar_vrt_desc];

        std::cout << "Edge:" << edge->getID() << std::endl;
        auto edge_Qres = edge->getLastQ(t_step - 1);
        auto edge_Pres = edge->getLastP(t_step - 1);

        cout << "\tLast Q: " << edge_Qres[0] << endl;
        cout << "\tLast P: " << edge_Pres[0] << endl;
    }

    return 0;
}