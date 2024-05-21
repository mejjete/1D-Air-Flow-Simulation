#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <string>
#include <type_traits>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <algorithm>
#include <mpi.h>
#include <omp.h>
#include "network_hybrid.hpp"

int mpi_size, mpi_rank;
int t_step, total_edges, total_vertices, dx;
double h, rho, a;

// Return VertexProperty object for each MPI process
VertexProperty readJsonNetwork(boost::property_tree::ptree &, int);

int main(int argc, char **argv)
{
    using namespace boost;
    using namespace std;

    MPI_Init(&argc, &argv);
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
    int edges = 0;

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
    }

    return vertex;
}