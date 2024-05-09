#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <string>
#include <type_traits>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <tuple>
#include <algorithm>
#include <mpi.h>
#include <omp.h>
#include "network_par.hpp"

// For argument message passing between MPI processes
typedef struct _mpi_message
{
    int ID;
    int s_step;
    int head;
    int tail;
    double alpha;
    double beta;
    double gamma;
    double h;
} _mpi_message_t;

static_assert(std::is_trivial<_mpi_message_t>::value,
        "_mpi_message_t must be trivially constructible");

int main(int argc, char **argv)
{
    using namespace boost;
    using namespace std;

    int mpi_size, mpi_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    int t_step;
    double h;
    EdgeProperty *edge;
    std::vector<VertexProperty> vertex_array;
    std::vector<_mpi_message_t> arguments;

    if(mpi_rank == 0)
    {
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
        h = json_network.get<double>("h");
        t_step = json_network.get<int>("t_step");
        const int total_vert = json_network.get<int>("vertices");
        const int total_edg = json_network.get<int>("edges");
        const double rho = json_network.get<double>("rho");
        const int a = json_network.get<int>("a");

        auto default_out = cout.flags();
        cout << "--------------------------------------------------" << endl;
        cout << "Simulation parameters" << endl;
        cout << "--------------------------------------------------" << endl;
        cout.flags(std::ios::left);
        cout << setw(12) << "dx:" << setw(10) << dx << endl;
        cout << setw(12) << "dt:" << setw(12) << h << endl;
        cout << setw(12) << "time:" << setw(12) << t_step << endl;
        cout << setw(12) << "vertices: " << setw(12) << total_vert << endl;
        cout << setw(12) << "edges: " << setw(12) << total_edg << endl;
        cout << setw(12) << "processes: " << mpi_size << setw(12) << endl;
        cout << "--------------------------------------------------" << endl;
        cout.flags(default_out);

        int total_edges = 0;

        if(mpi_size != total_edg)
        {
            std::cerr << "Total number of processes does not equal total number of edges... Exiting" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        /*                  EDGE INITIALIZATION              */
        boost::property_tree::ptree arrayEdges = json_network.get_child("edge");
        
        vector<double> vertex_gamma(total_vert, 0.0);
        vector<vector<int>> vertex_incoming(total_vert, vector<int>(0, 0));
        vector<vector<int>> vertex_outcoming(total_vert, vector<int>(0, 0));

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

            // Each vertex gets biggest gamma among incoming edges
            if(vertex_gamma[tail] < gamma)
                vertex_gamma[tail] = gamma;
            
            /**
             *  Each vertex contains a list of edge's IDs that represent incoming
             *  and outcoming edges. For a given vertex, outcoming edges is an edge
             *  that has given vertex as a head, and incoming edges is an edges that 
             *  has given vertex as a tail.
            */
            vertex_incoming[tail].push_back(id);
            vertex_outcoming[head].push_back(id);

            arguments.push_back({id, s_step, head, tail, alpha, beta, gamma, h});
            total_edges++;
        }
        
        if(total_edges != total_edg)
        {
            std::cerr << "Total number of edges does not match. File is corrupted\n";
            exit(EXIT_FAILURE);
        }

        /*                  VERTEX INITIALIZATION                */
        for(int i = 0; i < total_vert; i++)
            vertex_array.push_back(VertexProperty(i, h, 0.0, vertex_gamma[i], vertex_incoming[i], vertex_outcoming[i]));

        boost::property_tree::ptree array_pressure = json_network.get_child("vens");
        for(auto &child : array_pressure)
        {
            int vertex_ID = child.second.get<int>("vertex_id");
            double pressure = child.second.get<double>("pressure");
            vertex_array[vertex_ID].setDefault();
            vertex_array[vertex_ID].setP(pressure);
        }
    }

    // Broadcast t_step and h to every process
    MPI_Bcast(&t_step, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&h, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);  

    // Send edge's coefficients to a corresponding process
    if(mpi_rank == 0)
    {
        for(int i = 1; i < mpi_size; i++)
            MPI_Send(&arguments[i], sizeof(_mpi_message_t), MPI_CHAR, i, 0, MPI_COMM_WORLD);
        
        _mpi_message_t msg = arguments[0];
        edge = new EdgeProperty(msg.ID, msg.s_step, msg.head, msg.tail, msg.alpha, msg.beta, msg.gamma, msg.h);
    }
    else 
    {
        _mpi_message_t msg;
        MPI_Status status;
        MPI_Recv(&msg, sizeof(_mpi_message_t), MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
        edge = new EdgeProperty(msg.ID, msg.s_step, msg.head, msg.tail, msg.alpha, msg.beta, msg.gamma, msg.h);
    }


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
    std::vector<VertexDebug> vertex_debug;
    EdgeDebug edge_debug("E" + std::to_string(edge->getID()), edge);

    // Initialize vertex debug handlers
    for(auto &vert : vertex_array)
        vertex_debug.push_back(VertexDebug("V" + std::to_string(vert.getID()), &vert));

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

    int s_step = edge->getSteps();

    for(int i = 0; i < t_step; i++)
    {
        double source_P = 0.0;
        double target_P = 0.0;

        // Do boundary exchange
        if(mpi_rank == 0)
        {
            for(auto &vert : vertex_array)
            {
                double flow = 0.0;

                // Receive flows from both incoming and outcoming edges
                if(i != 0)
                {
                    double temp = 0.0;

                    for(auto ed : vert.getOutcoming())
                    {
                        MPI_Recv(&temp, 1, MPI_DOUBLE, ed, FLOW, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        flow += temp;
                    }

                    for(auto ed : vert.getIncoming())
                    {
                        MPI_Recv(&temp, 1, MPI_DOUBLE, ed, FLOW, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        flow += temp;
                    }
                }

                // Calculate flow in vertex
                double pressure = vert.calculateP(i, flow);

                // Send pressure to all incoming and outcoming edges
                for(auto ed : vert.getOutcoming())
                    MPI_Send(&pressure, 1, MPI_DOUBLE, ed, PRESSURE, MPI_COMM_WORLD);

                for(auto ed : vert.getIncoming())
                    MPI_Send(&pressure, 1, MPI_DOUBLE, ed, PRESSURE, MPI_COMM_WORLD);
                
                #ifdef DEBUG
                // Write down time points for each vertex separately
                if((i % (t_step / 3840)) == 0)
                {
                    auto &curr_vert_debug = vertex_finder(vert.getID());
                    curr_vert_debug.serialize(i * h);
                }
                #endif
            }
        }

        /**
         *  MPI guarantees to recieve messages in order they are sent.
         *  
         *  Additionally, we must guarantee that 2 messages has been sent in right order,
         *  so following calls to MPI_Recv do not yeild any delays.
         */
        MPI_Recv(&source_P, 1, MPI_DOUBLE, 0, PRESSURE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&target_P, 1, MPI_DOUBLE, 0, PRESSURE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);      

        edge->setSourceP(i, source_P);
        edge->setTargetP(i, target_P);

        /*************************************************************/

        // Loop over spatial steps for flow
        for(int k = 1; k < s_step; k++)
            edge->calculateQ(i, k);

        // Loop over spatial steps for pressure
        for(int k = 1; k < s_step - 1; k++)
            edge->calculateP(i, k);

        /*************************************************************/
        
        // Source flow goes to source vertex with inverse sign
        double source_Q = -edge->getLastQ(i, 1);

        // Target flow goes to target vertex with normal sign
        double target_Q = edge->getLastQ(i, s_step - 1);

        // Send incoming and outcoming flows to source and target vertices
        MPI_Send(&source_Q, 1, MPI_DOUBLE, 0, FLOW, MPI_COMM_WORLD);
        MPI_Send(&target_Q, 1, MPI_DOUBLE, 0, FLOW, MPI_COMM_WORLD);

        #ifdef DEBUG
        // Write down time points for each edge separately
        if((i % (t_step / 3840)) == 0)
            edge_debug.serialize(i * h);
        #endif
    }

    printf("Edge %d\n", edge->getID());
    printf("\tQ: %f\n", edge->getLastQ(t_step, s_step - 1));
    printf("\tP: %f\n", edge->getLastP(t_step, s_step - 1));

    delete edge;
    MPI_Finalize();
    return 0;
}