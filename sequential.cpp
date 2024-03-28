#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <string>
#include <omp.h>
#include <boost/circular_buffer.hpp>
#include <boost/container/vector.hpp>

int main() 
{
    using namespace std;

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
    const double alpha = 0.093;
    const double beta = 0.004;
    const double gamma = 728.6;

    /**
     * h      - discretization step in Euler method
     * tmax   - the number of approximation in Euler method
     * Simulation duration time = tmax * h
    */
    const int tmax = 1800000;
    const double h = 0.0001;

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

#ifdef DEBUG     
    const int t_step = tmax;
    const int s_step = kmax;
#else
    const int t_step = 2;
    const int s_step = kmax;
#endif

    auto Q = new double[t_step][s_step] {{0.0}};
    auto P = new double[t_step][s_step] {{0.0}};

    // Set initial condition for pressure
    for(int i = 0; i < t_step; i++)
        P[i][s_step - 1] = -202.0;
    
    /**
     * In fact, we have two steps. The first step represents the approximation 
     * criteria in Euler's method (i) and the second one represents the actual
     * spatial step along the edge (k).
    */
    for (int i = 0; i < tmax; i++) 
    {
        int i_next = (i + 1) % t_step;
        int i_curr = i % t_step;

        // Loop over spatial steps for flow
        for (int k = 1; k < s_step; k++)
            Q[i_next][k] = Q[i_curr][k] + h * (alpha * (P[i_curr][k - 1] - P[i_curr][k]) - beta * Q[i_curr][k] * abs(Q[i_curr][k]));

        // Loop over steps for pressure, DO NOT calculate pressure for the last element as it is given as a boundary condition
        for (int k = 1; k < s_step - 1; k++)
            P[i_next][k] = P[i_curr][k] + h * (gamma * (Q[i_curr][k] - Q[i_curr][k + 1]));  
    }

    cout << "First 10 spatial steps of last temporal step\n";
    cout << "[1 - 10][" << tmax - 1 << "]";
    for(int i = 0; i < 10; i++)
        printf("%10.3f", Q[t_step - 1][i + 1]); 
    printf("\n");

#ifdef DEBUG
    // Generate text files for gnuplot
    FILE *Q_plot_s = fopen("Q_plots_s.txt", "w+");
    FILE *Q_plot_m = fopen("Q_plots_m.txt", "w+");
    FILE *Q_plot_e = fopen("Q_plots_e.txt", "w+");

    FILE *P_plot_s = fopen("P_plots_s.txt", "w+");
    FILE *P_plot_m = fopen("P_plots_m.txt", "w+");
    FILE *P_plot_e = fopen("P_plots_e.txt", "w+");

    for (int i = 1; i < t_step; i += 10)
    {
        fprintf(Q_plot_s, "%d %f\n", int(i * h), Q[i][1]);
        fprintf(P_plot_s, "%d %f\n", int(i * h), P[i][1]);

        fprintf(Q_plot_m, "%d %f\n", int(i * h), Q[i][s_step / 2]);
        fprintf(P_plot_m, "%d %f\n", int(i * h), P[i][s_step / 2]);  

        fprintf(Q_plot_e, "%d %f\n", int(i * h), Q[i][s_step - 1]);
        fprintf(P_plot_e, "%d %f\n", int(i * h), P[i][s_step - 1]);
    }

    fclose(Q_plot_s);
    fclose(Q_plot_m);
    fclose(Q_plot_e);
    fclose(P_plot_s);
    fclose(P_plot_m);
    fclose(P_plot_e);
#endif
    return 0;
}
