#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>

using namespace std;

int main() 
{
    const double F = 4.67;
    const double r = 0.00105;
    const double rho = 1.255;
    const double air = 330.22;

    const int l = 1925;               // Edge length
    const int dx = 40;                // Spatial step
    const int kmax = l / dx;          // Number of spatial steps

    /**
     * Coefficients alpha beta and gamma are different for different 
     * edges. It depends on the edge's lenght, cross-sectional area etc.
    */
    const double alpha = F / (rho * dx);
    const double beta = (F * r) / rho;
    const double gamma = (rho * air * air) / (F * dx);

    /**
     * h      - discretization step in Euler method
     * tmax   - the number of approximation in Euler method
     * Changing of any of them might result in NaN error.
    */
    const int tmax = 100000;
    const double h = 0.001;

    auto default_out = cout.flags();
    cout << "--------------------------------------------------" << endl;
    cout << "Simulation parameters" << endl;
    cout << "--------------------------------------------------" << endl;
    cout.flags(std::ios::left);
    cout << setw(10) << "l:" << setw(10) << l << setw(5) << "|";
    cout << setw(10) << "alpha: " << alpha << setw(10) << endl;
    cout << setw(10) << "dx:" << setw(10) << dx << setw(5) << "|";
    cout << setw(10) << "beta: " << beta << setw(10) << endl;
    cout << setw(10) << "kmax:" << setw(10) << std::setw(10) << kmax << setw(5) << "|";
    cout << setw(10) << "gamma: " << gamma << setw(10) << endl;
    cout << setw(10) << "tmax:" << setw(10) << std::setw(10) << tmax << setw(5) << "|" << endl;
    cout << setw(10) << "h:" << setw(10) << std::setw(10) << h << setw(5) << "|" << endl;
    cout << "--------------------------------------------------" << endl;
    cout.flags(default_out);

    // Stores the approximated values for each step 
    auto Q = new double[tmax][kmax];
    auto P = new double[tmax][kmax];

    // Prepare initial condition
    for (int i = 0; i < tmax; i++) 
    {
        // Loop over spatial steps
        for (int k = 0; k < kmax; k++) 
        {
            // Boundary condition for the last spatial step
            if (k == kmax - 1) 
            {
                Q[i][k] = 0;
                P[i][k] = -202.0;   // P = -202 when Q = 10
            } else 
            {
                Q[i][k] = 0;
                P[i][k] = 0;
            }
        }
    }

    /**
     * In fact, we have two steps. The first step represents the approximation 
     * criteria in Euler's method (tmax) and the second one represents the actual
     * spatial step along the edge (dx).
    */
    for (int i = 2; i < tmax; i++) 
    {
        // Loop over spatial steps for flow
        for (int k = 1; k < kmax; k++)
            Q[i + 1][k] = Q[i][k] + h * (alpha * (P[i][k - 1] - P[i][k]) - beta * Q[i][k] * abs(Q[i][k]));

        // Loop over steps for pressure, DO NOT calculate pressure for the last element as it is given as a boundary condition
        for (int k = 1; k < kmax - 1; k++)
            P[i + 1][k] = P[i - 2][k] + h * (gamma * (Q[i][k] - Q[i][k + 1]));
    }

    cout << "First 10 spatial steps of last temporal step\n";
    cout << "[1 - 10][" << tmax - 1 << "]";
    for(int i = 0; i < 10; i++)
        printf("%10.3f", Q[tmax - 1][i + 1]); 

    return 0;
}
