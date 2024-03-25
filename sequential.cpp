#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <string>
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

    // Just a thumbs for future use
    const int t_step = tmax;
    const int s_step = kmax;

    auto Q = new double[t_step][s_step] {0.0};
    auto P = new double[t_step][s_step] {0.0};

    // Helper function to get k depending on the Q and P context
    auto ind = [](int offset = 0, int k) { return k - offset; };

    // Set initial condition for pressure
    for(int i = 0; i < t_step; i++)
        P[i][s_step - 1] = -202.0;
    
    /**
     * In fact, we have two steps. The first step represents the approximation 
     * criteria in Euler's method (i) and the second one represents the actual
     * spatial step along the edge (k).
    */
    for (int i = 0; i < t_step; i++) 
    {
        // Loop over spatial steps for flow
        for (int k = 1; k < s_step; k++)
            Q[i + 1][k] = Q[i][k] + h * (alpha * (P[i][k - 1] - P[i][k]) - beta * Q[i][k] * abs(Q[i][k]));

        // Loop over steps for pressure, DO NOT calculate pressure for the last element as it is given as a boundary condition
        for (int k = 1; k < s_step - 1; k++)
            P[i + 1][k] = P[i][k] + h * (gamma * (Q[i][k] - Q[i][k + 1]));
    }

    cout << "First 10 spatial steps of last temporal step\n";
    cout << "[1 - 10][" << tmax - 1 << "]";
    for(int i = 0; i < 10; i++)
        printf("%10.3f", Q[t_step - 1][i + 1]); 
    printf("\n");

#ifdef DEBUG    
    FILE *gnuplot = popen("gnuplot -persist", "w");
    if (gnuplot == NULL) 
    {
        perror("Error opening pipe to GNU Plot");
        return -1;
    }

    fprintf(gnuplot, "set terminal png size 1920,1080\n");
    fprintf(gnuplot, "set output 'Q.png'\n");

    fprintf(gnuplot, "set xlabel 'Time (seconds)'\n");
    fprintf(gnuplot, "set ylabel 'Q (flow)'\n");

    std::string plot_option = "plot '-' with lines title ";
    plot_option += "'Step: " + std::to_string(s_step - 1) + "' linecolor rgb 'red',";
    plot_option += "'-' with line title 'Step: " + std::to_string(s_step / 2) + "' linecolor rgb 'blue',"; 
    plot_option += "'-' with line title 'Step: " + std::to_string(1) + "' linecolor rgb 'green'\n"; 
    fprintf(gnuplot, plot_option.c_str());

    // Output last step
    for (int i = 1; i < t_step; i += 10) 
        fprintf(gnuplot, "%d %f\n", int(i * h), Q[i][s_step - 1]);
    fprintf(gnuplot, "e\n");

    // Output middle step
    for (int i = 1; i < t_step; i += 10) 
        fprintf(gnuplot, "%d %f\n", int(i * h), Q[i][s_step / 2]);
    fprintf(gnuplot, "e\n");

    // Output first step
    for (int i = 1; i < t_step; i += 10) 
        fprintf(gnuplot, "%d %f\n", int(i * h), Q[i][1]);
    fprintf(gnuplot, "e\n");

    if (pclose(gnuplot) == -1) 
        perror("Error closing pipe to GNU Plot");
#endif

    return 0;
}
