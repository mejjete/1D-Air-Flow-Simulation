#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>

using namespace std;

int main() 
{
    int l = 1925;             // Edge length
    int dx = 40;               // Spatial step
    int dt = 100000;          // Euler's step
    int kmax = l / dx;        // Number of spatial steps

    // Parameters
    double alpha = 0.093;
    double beta = 0.004;
    double gamma = 728.6;
    double h = 0.007;

    // Arrays to store data
    double Q[kmax][dt];
    double P[kmax][dt];

    // Loop over time steps
    for (int i = 0; i < dt; i++) 
    {
        // Loop over spatial steps
        for (int k = 0; k < kmax; k++) 
        {
            // Boundary condition for the last spatial step
            if (k == kmax - 1) 
            {
                Q[k][i] = 0;
                P[k][i] = -202.0;   // P = -202 when Q = 10
            } else 
            {
                Q[k][i] = 0;
                P[k][i] = 0;
            }
        }
    }

    // Main calculation loop
    for (int i = 2; i < dt; i++) 
    {
        // Loop over spatial steps
        for (int k = 1; k < kmax; k++) 
        {
            Q[k][i + 1] = Q[k][i] + h * (alpha * (P[k - 1][i] - P[k][i]) - beta * Q[k][i] * abs(Q[k][i]));

            // DO NOT calculate pressure for the last element such as it is given by default
            if (k != kmax - 1) 
                P[k][i + 1] = P[k][i - 2] + h * (gamma * (Q[k][i] - Q[k + 1][i]));
        }
    }
    
    // Print out the first 10 of Q for each spatial step
    std::cout << "FIRST 10 VALUES\n";
    for(int i = 0; i < kmax; i++)
    {
        std::cout.width(5);
        std::cout.flags(std::ios::left);
        std::cout << i << ": ";
   
        for(int k = 0; k < 10; k++)
            std::cout << std::setw(10) << Q[i][k];
        std::cout << std::endl;
    }

    // Print out the first 10 of Q for each spatial step
    std::cout << "LAST 10 VALUES:\n";
    for(int i = 1; i < kmax; i++)
    {
        std::cout.width(5);
        std::cout.flags(std::ios::left);
        std::cout << i << ": ";

        for(int k = 0; k < 10; k++)
            std::cout << std::setw(10) << Q[i][dt - 10 + k];
        std::cout << std::endl;
    }


    return 0;
}
