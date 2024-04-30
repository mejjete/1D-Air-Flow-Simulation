/**
 *  Parallel implementation using only OpenMP.
 * 
 *  There's no actual loop scheduling. We have two hard-coded edges
 *  that is executed by 2 threads in parallel using sections directive.
 *  It just compares two results between single-threaded and multi-threaded
 *  instances to check any difference. 
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

    EdgeProperty *e3 = new EdgeProperty(3, 49, t_step, 0.093, 0.004, 728.6, 0.0001);
    e3->setSourceP(0, 0);
    e3->setTargetP(0, -202);
    e3->setTargetP(1, -202);

    EdgeProperty *e4 = new EdgeProperty(4, 35, t_step, 0.12, 0.009, 569.1, 0.0001);
    e4->setSourceP(0, 0);
    e4->setTargetP(0, -594.76);
    e4->setTargetP(1, -594.76);

    EdgeProperty *network[] = {e3, e4};
    int elements = sizeof(network) / sizeof(EdgeProperty *);

    double start = omp_get_wtime();
    #pragma omp parallel
    {
        #pragma omp single
        std::cout << "Threads: " << omp_get_num_threads() << std::endl;

        for(int i = 0; i < t_step; i++)
        {
            #pragma omp sections
            {
                #pragma omp section
                {
                    calculateQ(e3, i);
                    calculateP(e3, i);
                }

                #pragma omp section
                {
                    calculateQ(e4, i);
                    calculateP(e4, i);
                }
            }
        }
    }
    double end = omp_get_wtime();
    
    std::cout << "E3-Q: " << e3->getLastQ(t_step)[0] << std::endl;
    std::cout << "E3-P: " << e3->getLastP(t_step)[0] << std::endl;

    std::cout << "E4-Q: " << e4->getLastQ(t_step)[0] << std::endl;
    std::cout << "E4-P: " << e4->getLastP(t_step)[0] << std::endl;
    std::cout << std::endl << "Time: " << end - start << std::endl;
}
