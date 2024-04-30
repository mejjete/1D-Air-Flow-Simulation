/**
 *  Parallel implementation using hybrid scheme MPI + OpenMP.
 * 
 *  MPI holds and computes each edge, and OpenMP with 2 threads 
 *  computes Q (flow) and P(pressure) in parallel using sections.
 *  
 *  Entire network is hard-coded and does not exchange boundaries as it
 *  should. Instead, it just computes each edge separately.
*/

#include <iostream>
#include <cstddef>
#include <vector>
#include <omp.h>
#include <mpi.h>
#include "../network.hpp"

using namespace std;

void calculateQ(EdgeProperty *edge, int t_step)
{
    for(int k = 1; k < edge->getSteps(); k++)
        edge->calculateQ(t_step, k);
}

void calculateP(EdgeProperty *edge, int t_step)
{
    for(int k = 1; k < edge->getSteps() - 1; k++)
        edge->calculateP(t_step, k);
}

int main(int argc, char **argv)
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if(size != 8)
    {
        std::cerr << "Invalid number of processes... Exiting\n";
        abort();
    }

    const int t_step = 180000;
    EdgeProperty *private_edge;

    // Each process should allocate it's memory region by itself
    if(rank == 0)
    {
        private_edge = new EdgeProperty(1, 20, t_step, 0.13, 0.00252, 523.6, 0.0001);
        private_edge->setTargetP(0, -1985.0);
        private_edge->setTargetP(1, -1985.0);
    }
    else if(rank == 1)
    {
        private_edge = new EdgeProperty(2, 4, t_step, 0.13, 0.00252, 523.6, 0.0001);
        private_edge->setTargetP(0, -60.63);
        private_edge->setTargetP(1, -60.63);
    }
    else if(rank == 2)
    {
        private_edge = new EdgeProperty(3, 49, t_step, 0.093, 0.004, 728.6, 0.0001);
        private_edge->setTargetP(0, -202);
        private_edge->setTargetP(1, -202);
    }
    else if(rank == 3)
    {
        private_edge = new EdgeProperty(4, 35, t_step, 0.12, 0.009, 569.1, 0.0001);
        private_edge->setTargetP(0, -594.76);
        private_edge->setTargetP(1, -594.76);
    }
    else if(rank == 4)
    {
        private_edge = new EdgeProperty(5, 43, t_step, 0.09, 0.004, 756.2, 0.0001);
        private_edge->setTargetP(0, -710.14);
        private_edge->setTargetP(1, -710.14);
    }
    else if(rank == 5)
    {
        private_edge = new EdgeProperty(6, 24, t_step, 0.1, 0.0055, 684.7, 0.0001);
        private_edge->setTargetP(0, -718.72);
        private_edge->setTargetP(1, -718.72);
    }
    else if(rank == 6)
    {
        private_edge = new EdgeProperty(7, 20, t_step, 0.13, 0.00252, 523.6, 0.0001);
        private_edge->setTargetP(0, -60.63);
        private_edge->setTargetP(1, -60.63);
    }
    else if(rank == 7)
    {
        private_edge = new EdgeProperty(8, 17, t_step, 0.13, 0.00252, 523.6, 0.0001);
        private_edge->setTargetP(0, -1599.4);
        private_edge->setTargetP(1, -1599.4);
    }

    double start = omp_get_wtime();
    #pragma omp parallel firstprivate(size, rank)
    {
        omp_set_num_threads(2);

        #pragma omp single
        {
            if(rank == 0)
                std::cout << "Threads: " << omp_get_num_threads() << std::endl;
        }

        for(int i = 0; i < t_step; i++)
        {
            #pragma omp sections
            {
                #pragma omp section
                    calculateQ(private_edge, i);

                #pragma omp section
                    calculateP(private_edge, i);
            }
        }

        // Boundary exchange
        #pragma omp barrier
        double bnd = omp_get_wtime() - start;
    }
    double end = omp_get_wtime();
    
    std::cout << "Edge_" << private_edge->getID() << std::endl;
    std::cout << "\tQ: " << private_edge->getLastQ(t_step)[0] << std::endl;
    std::cout << "\tP: " << private_edge->getLastP(t_step)[0] << std::endl;

    MPI_Finalize();
    return 0;
}
