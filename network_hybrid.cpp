#include <cstdio>
#include <iostream>
#include <iomanip>
#include <string>
#include <set>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <mpi.h>
#include <omp.h>
#include <vector>
#include "network_hybrid.hpp"

int mpi_size, mpi_rank;
int t_step, total_edges, total_vertices, dx;
double h, rho, a;

// Return VertexProperty object in each MPI process
VertexProperty readJsonNetwork(boost::property_tree::ptree &, int);

int main(int argc, char **argv)
{
    using namespace boost;
    using namespace std;

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    
    if(provided < MPI_THREAD_MULTIPLE) 
    {
        std::cerr << "MPI does not provide required thread support" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    // User must specify JSON configuration file for entire network
    if(argc <= 1)
        throw std::runtime_error("No input file exist\n");

    boost::property_tree::ptree json_network;
    boost::property_tree::read_json(argv[1], json_network);
    const std::string network_name = json_network.get<std::string>("name");
    dx = json_network.get<int>("dx");
    h = json_network.get<double>("h");
    t_step = json_network.get<int>("t_step");
    total_vertices = json_network.get<int>("vertices");
    total_edges = json_network.get<int>("edges");
    rho = json_network.get<double>("rho");
    a = json_network.get<int>("a");

    if(mpi_size != total_vertices)
        throw std::runtime_error("Number of processes must be equal to number of vertices\n"); 

    if(mpi_rank == 0)
    {
        auto default_out = cout.flags();
        cout << "--------------------------------------------------" << endl;
        cout << "Simulation parameters" << endl;
        cout << "--------------------------------------------------" << endl;
        cout.flags(std::ios::left);
        cout << setw(12) << "dx:" << setw(10) << dx << endl;
        cout << setw(12) << "dt:" << setw(12) << h << endl;
        cout << setw(12) << "time:" << setw(12) << t_step << endl;
        cout << setw(12) << "vertices: " << setw(12) << total_vertices << endl;
        cout << setw(12) << "edges: " << setw(12) << total_edges << endl;
        cout << setw(12) << "processes: " << mpi_size << setw(12) << endl;
        cout << "--------------------------------------------------" << endl;
        cout.flags(default_out);
    }

    VertexProperty vertex = readJsonNetwork(json_network, mpi_rank);

    #if 1
    // Team size is equal to number of outcoming edges in dedicated vertex
    int team_size = vertex.getEdgeCount();
    omp_set_num_threads(team_size);

    #pragma omp parallel
    {
        #pragma omp single
        {
            if(omp_get_num_threads() != team_size && team_size != 0)
                throw std::runtime_error("Number of threads != number of edges associated with vertex\n");
        }

        // Head vertex has no incoming edges
        if(team_size > 0 && vertex.getID() != 0)
        {
            int thread_id = omp_get_thread_num();
                EdgeProperty &edge = vertex.getEdge(thread_id);
            int s_step = edge.getSteps();

            // Main calculation
            for(int i = 0; i < t_step; i++)
            {
                for(int k = 1; k < s_step; k++)
                    edge.calculateQ(i, k);
                
                for(int k = 1; k < s_step - 1; k++)
                    edge.calculateP(i, k);

                #pragma omp master
                {
                    // Exchange parameters
                    // double source_Q;
                }

                // Wait till master thread completes the communication so we can continue calculation
                #pragma omp barrier
            }

            #pragma omp critical
            {
                printf("Edge %d\n", edge.getID());
                printf("\tQ: %f\n", edge.getLastQ(t_step, s_step - 1));
                printf("\tP: %f\n", edge.getLastP(t_step, s_step - 1));
            }
        }
    }

    #endif

    // printf("Vertex: %d\n", vertex.getID());
    // for(auto iter : vertex.getEdgeGroups())
    // {
    //     std::cout << "\tGroup: ";
    //     for(auto edge : iter.getEdges())
    //         std::cout << edge.getID() << " ";
    //     std::cout << std::endl;
    // } 

    // printf("Vertex: %d\n", vertex.getID());
    // printf("\tExpecting messages: %d\n", vertex.getVertexNum());

    MPI_Finalize();
    return 0;
}

VertexProperty readJsonNetwork(boost::property_tree::ptree &json_network, int rank)
{
    /*                  VERTEX INITIALIZATION                */
    VertexProperty vertex(rank, h, 0, 0);

    // Set default pressure for each vertex that has vens on it
    boost::property_tree::ptree array_pressure = json_network.get_child("vens");
    for(auto &child : array_pressure)
    {
        int vertex_ID = child.second.get<int>("vertex_id");
        double pressure = child.second.get<double>("pressure");

        if(vertex_ID == rank)
        {
            vertex.setDefault();
            vertex.setP(pressure);
        }
    }

    /*                  EDGE INITIALIZATION              */
    boost::property_tree::ptree arrayEdges = json_network.get_child("edge");
    std::set<std::pair<int, int>> out_edges;
    std::set<int> vert;

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

        if(tail == rank)
        {
            EdgeProperty edge(id, s_step, head, tail, alpha, beta, gamma, h);
            vertex.addEdge(edge);
        }

        /**
         *  Parallel Edge Reduction mechanism requires us to know hom much messages we are gonna
         *  receive. Because some outcoming edges can be grouped together and sent as a 
         *  1 single message, we have to make sure we receive all messages by tracking message amount.
        */
        if(head == rank)
        {
            out_edges.insert(std::make_pair(head, tail));
            vert.insert(tail);
        }
    }

    // Add adjacent vertices
    for(auto _[[maybe_unused]] : vert)
        vertex.addVertex();

    vertex.setOutMsg(out_edges.size());
    return vertex;
}

void VertexProperty::addEdge(EdgeProperty edge)
{
    // For a given vertex we take the biggest gamma among its edges
    if(edge.getGamma() > gamma)
        gamma = edge.getGamma();

    edge.setTargetP(0, getP(0));
    edge.setTargetP(1, getP(0));

    int edge_head = edge.getHead();

    auto iter = std::find_if(edge_groups.begin(), edge_groups.end(), [edge_head](EdgeGroup &grp) 
        { return edge_head == grp.getVertex(); });
    
    if(iter == edge_groups.end())
    {
        EdgeGroup new_group(edge_head);
        new_group.addEdge(edge);
        edge_groups.push_back(new_group);
    }
    else
        iter->addEdge(edge);
}

EdgeProperty &VertexProperty::getEdge(size_t edge_ID)
{
    for(auto &group : edge_groups)
    {
        if(edge_ID >= group.getEdges().size())
        {
            edge_ID -= group.getEdges().size();
            continue;
        }
        else 
            return group.getEdges()[edge_ID];
    }

    throw std::runtime_error("Edge: " + std::to_string(edge_ID) + " does not exist in vertex " 
        + std::to_string(getID()));
}

int VertexProperty::getEdgeCount()
{
    int count = 0;

    for(auto &group : edge_groups)
        count += group.getEdges().size();
    
    return count;
}
