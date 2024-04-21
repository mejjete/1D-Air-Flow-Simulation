#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <string>
#include <omp.h>
#include <boost/circular_buffer.hpp>
#include <boost/container/vector.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

int main(int argc, char **argv) 
{
    using namespace std;

    if(argc <= 1)
    {
        std::cerr << "No input file exist\n";
        exit(EXIT_FAILURE);
    }

    boost::property_tree::ptree edge;
    boost::property_tree::read_json(argv[1], edge);

    const std::string edge_name = edge.get<std::string>("name");
    const int length = edge.get<int>("length");
    const int dx = edge.get<int>("dx");
    const double F = edge.get<double>("F");
    const double r = edge.get<double>("r");
    const double rho = edge.get<double>("rho");
    const int a = edge.get<int>("a");
    const double alpha = F / (rho * dx);
    const double beta = (F * r) / rho; 
    const double gamma = (rho * a * a) / (F * dx);
    const double h = edge.get<double>("h");
    const double P_init = edge.get<double>("P_init");
    const int t_step = edge.get<int>("t_step");
    const int s_step = (int)(ceil((double)length / (double)dx));

    auto default_out = cout.flags();
    cout << "--------------------------------------------------" << endl;
    cout << "Simulation parameters" << endl;
    cout << "--------------------------------------------------" << endl;
    cout.flags(std::ios::left);
    cout << setw(10) << "l:" << setw(10) << length << setw(5) << "|";
    cout << setw(10) << "alpha: " << alpha << setw(10) << endl;
    cout << setw(10) << "dx:" << setw(10) << dx << setw(5) << "|";
    cout << setw(10) << "beta: " << beta << setw(10) << endl;
    cout << setw(10) << "M:" << setw(10) << s_step << setw(5) << "|";
    cout << setw(10) << "gamma: " << gamma << setw(10) << endl;
    cout << setw(10) << "time:" << setw(10) << t_step << setw(5) << "|";
    cout << setw(10) << "P init: " << P_init << setw(10) << endl;
    cout << setw(10) << "dt:" << setw(10) << h << setw(5) << "|" << endl;
    cout << "--------------------------------------------------" << endl;
    cout.flags(default_out);

#ifdef DEBUG     
    // Generate text files for gnuplot
    FILE *Q_plot_s = fopen((edge_name + "_Qplots_s.txt").c_str(), "w+");
    FILE *Q_plot_m = fopen((edge_name + "_Qplots_m.txt").c_str(), "w+");
    FILE *Q_plot_e = fopen((edge_name + "_Qplots_e.txt").c_str(), "w+");

    FILE *P_plot_s = fopen((edge_name + "_Pplots_s.txt").c_str(), "w+");
    FILE *P_plot_m = fopen((edge_name + "_Pplots_m.txt").c_str(), "w+");
    FILE *P_plot_e = fopen((edge_name + "_Pplots_e.txt").c_str(), "w+");
#endif

    /**
     *  Q - air flow; P - pressure
     *  They represent circular buffer where rows are organized in a circular maner
    */
    double **Q = new double*[2];
    double **P = new double*[2]; 

    for(int i = 0; i < 2; i++)
    {
        Q[i] = new double[s_step];
        P[i] = new double[s_step];

        // Set initial condition for pressure
        P[i][s_step - 1] = P_init;
    }
    
    /**
     * In fact, we have two steps. The first step represents the approximation 
     * criteria in Euler's method (i) and the second one represents the actual
     * spatial step along the edge (k).
    */
    for (int i = 0; i < t_step; i++) 
    {
        // As soon as we have circular buffer, we need to adjust pointers at each iteration
        int i_next = (i + 1) % 2;
        int i_curr = i % 2;

        // Loop over spatial steps for flow
        for (int k = 1; k < s_step; k++)
        {
            Q[i_next][k] = Q[i_curr][k] + h * (alpha * (P[i_curr][k - 1] - P[i_curr][k]) - beta * Q[i_curr][k] * abs(Q[i_curr][k]));

            #ifdef DEBUG
            /**
             *  Let's assume we are working on a 4k monitor with resolution 3840 x 2160. 
             *  Graphing 1 point per 1 pixel is enough.
            */
            if((i % (t_step / 3840)) == 0)
            {
                if(k == 1)
                    fprintf(Q_plot_s, "%.3f %f\n", i * h, Q[i_next][k]);
                else if(k == s_step / 2)
                    fprintf(Q_plot_m, "%.3f %f\n", i * h, Q[i_next][k]);
                else if(k == s_step - 1)
                    fprintf(Q_plot_e, "%.3f %f\n", i * h, Q[i_next][k]);
            }
            #endif 
        }

        // Loop over steps for pressure, DO NOT calculate pressure for the last element as it is given as a boundary condition
        for (int k = 1; k < s_step - 1; k++)
        {
            P[i_next][k] = P[i_curr][k] + h * (gamma * (Q[i_curr][k] - Q[i_curr][k + 1]));
            
            #ifdef DEBUG
            // For an explanation look at the previous DEBUG section
            if((i % (t_step / 3840)) == 0)
            {
                if(k == 1)
                    fprintf(P_plot_s, "%.3f %f\n", i * h, P[i_next][k]);
                else if(k == s_step / 2)
                    fprintf(P_plot_m, "%.3f %f\n", i * h, P[i_next][k]);  
                else if(k == s_step - 2)
                    fprintf(P_plot_e, "%.3f %f\n", i * h, P[i_next][k]);
            }
            #endif 
        }
    }

    printf("Q[%d][%d]: %.3f\n", t_step - 1, s_step - 1, Q[1][s_step - 1]);
    printf("P[%d][%d]: %.3f\n", t_step - 1, s_step - 2, P[1][s_step - 2]);

#ifdef DEBUG
    fclose(Q_plot_s);
    fclose(Q_plot_m);
    fclose(Q_plot_e);
    fclose(P_plot_s);
    fclose(P_plot_m);
    fclose(P_plot_e);
#endif
    
    for(int i = 0; i < 2; i++)
    {
        delete[] Q[i];
        delete[] P[i];
    }

    return 0;
}
