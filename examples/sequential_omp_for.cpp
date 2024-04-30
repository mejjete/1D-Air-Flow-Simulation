/**
 *  Parallel implementation using only OpenMP.
 * 
 *  In this case we have all test network hard-coded not only 2 edges.
 *  I used worksharing for directive to split the job.
*/

#include <iostream>
#include <cstddef>
#include <vector>
#include <omp.h>
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

int main()
{
    const int t_step = 3600000;

    EdgeProperty *e1 = new EdgeProperty(1, 20, t_step, 0.13, 0.00252, 523.6, 0.0001);
    e1->setTargetP(0, -1985.0);
    e1->setTargetP(1, -1985.0);

    EdgeProperty *e2 = new EdgeProperty(2, 4, t_step, 0.13, 0.00252, 523.6, 0.0001);
    e2->setTargetP(0, -60.63);
    e2->setTargetP(1, -60.63);

    EdgeProperty *e3 = new EdgeProperty(3, 49, t_step, 0.093, 0.004, 728.6, 0.0001);
    e3->setTargetP(0, -202);
    e3->setTargetP(1, -202);

    EdgeProperty *e4 = new EdgeProperty(4, 35, t_step, 0.12, 0.009, 569.1, 0.0001);
    e4->setTargetP(0, -594.76);
    e4->setTargetP(1, -594.76);

    EdgeProperty *e5 = new EdgeProperty(5, 43, t_step, 0.09, 0.004, 756.2, 0.0001);
    e5->setTargetP(0, -710.14);
    e5->setTargetP(1, -710.14);

    EdgeProperty *e6 = new EdgeProperty(6, 24, t_step, 0.1, 0.0055, 684.7, 0.0001);
    e6->setTargetP(0, -718.72);
    e6->setTargetP(1, -718.72);

    EdgeProperty *e7 = new EdgeProperty(7, 20, t_step, 0.13, 0.00252, 523.6, 0.0001);
    e7->setTargetP(0, -60.63);
    e7->setTargetP(1, -60.63);

    EdgeProperty *e8 = new EdgeProperty(8, 17, t_step, 0.13, 0.00252, 523.6, 0.0001);
    e8->setTargetP(0, -1599.4);
    e8->setTargetP(1, -1599.4);

    // EdgeProperty *network[] = {e1, e2, e3, e4, e5, e6, e7, e8};
    EdgeProperty *network[] = {e3};
    int elements = sizeof(network) / sizeof(EdgeProperty *);

    double start = omp_get_wtime();
    #pragma omp parallel
    {
        #pragma omp single
        std::cout << "Threads: " << omp_get_num_threads() << std::endl;

        for(int i = 0; i < t_step; i++)
        {
            #pragma omp for
            for(int e = 0; e < elements; e++)
            {
                calculateQ(network[e], i);
                calculateP(network[e], i);
            }
        }

        // Boundary exchange
        #pragma omp barrier
        double bnd = omp_get_wtime() - start;
    }
    double end = omp_get_wtime();
    
    for(int i = 0; i < elements; i++)
    {
        std::cout << "Edge_" << network[i]->getID() << std::endl;
        std::cout << "\tQ: " << network[i]->getLastQ(t_step)[0] << std::endl;
        std::cout << "\tP: " << network[i]->getLastP(t_step)[0] << std::endl;
    }

    std::cout << std::endl << "Time: " << end - start << std::endl;
}
