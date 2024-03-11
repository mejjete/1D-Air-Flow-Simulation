#include <iostream>
#include <vector>
#include <iomanip>
#include <math.h>

using namespace std;

int main() 
{
    // Length
    int l = 1925;
    // Spatial step
    int dx = 40;
    // Temporal step
    int dt = 2500;
    // Number of spatial steps
    int kmax = l / dx;

    // Initializing parameters
    vector<double> list(dt);

    // Arrays to store data
    double Q[kmax][dt];
    double P[kmax][dt];

    int i;

    // Parameters
    double alpha = 0.093;
    double beta = 0.004;
    double gamma = 728.6;
    double h = 0.05;
    double sum_of_h = 0;

    // Loop over time steps
    for (i = 0; i < dt; i++) 
    {
        // Loop over spatial steps
        for (int k = 0; k < kmax; k++) 
        {
            // Boundary condition for the last spatial step
            if (k == kmax - 1) 
            {
                Q[k][i] = 0;
                P[k][i] = -202.0; // P = -202 when Q = 10
            } else 
            {
                Q[k][i] = 0;
                P[k][i] = 0;
            }
        }
    }

    // Main calculation loop
    for (i = 2; i < dt - 2; i++) 
    {
        // Add current time step to the list
        sum_of_h = sum_of_h + h;
        list.push_back(sum_of_h);

		std::cout << "Q at step " << i;
		
        // Loop over spatial steps
        for (int k = 1; k < kmax; k++) 
        {
            // Boundary condition for the last spatial step
            if (k == kmax - 1) 
            {
                Q[k][i + 1] = Q[k][i] + h * (alpha * (P[k - 1][i] - P[k][i]) - beta * Q[k][i] * abs(Q[k][i]));
            } else 
            {
                // Calculate Q and P for non-boundary spatial steps
                Q[k][i + 1] = Q[k][i] + h * (alpha * (P[k - 1][i] - P[k][i]) - beta * Q[k][i] * abs(Q[k][i]));
                P[k][i + 1] = P[k][i - 2] + h * (gamma * (Q[k][i] - Q[k + 1][i]));
            }
        }
    }
    
    // Print out the first 10 values of the timestemp
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

    // Print out the last 10 values of the timestemp
    std::cout << "LAST 10 VALUES\n";
    for(int i = 0; i < kmax; i++)
    {
        std::cout.width(5);
        std::cout.flags(std::ios::left);
        std::cout << i << ": ";
 
        for(int k = 0; k < 10; k++)
            std::cout << std::setw(10) << Q[i][2489 + k];
        std::cout << std::endl;
    }

    return 0;
}
