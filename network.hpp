#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

// Property map for edge
class EdgeProperty
{
private:
    int length;
    int steps;
    int id;
    double alpha;
    double beta;
    double gamma;
    double h;

    // Circular buffers for flow and pressure
    std::vector<std::vector<double>> Q;
    std::vector<std::vector<double>> P;

public:
    // Edge initialization. The last parameter specifies the boundary condition for pressure
    EdgeProperty(int ID, int l, int s, int t_steps, double al, double bet, double gam, double hh)
        : id(ID), length(l), steps(s), alpha(al), beta(bet), gamma(gam), h(hh), 
        Q(2, std::vector<double>(t_steps, 0.0)), P(2, std::vector<double>(t_steps, 0.0))
    {};

    // One-spatial-step calculation of air flow
    double calculateQ(int t_step, int s_step)
    {
        int i_next = (t_step + 1) % 2;
        int i_curr = t_step % 2;
        int k = s_step;

        return Q[i_next][k] = Q[i_curr][k] + h * (alpha * (P[i_curr][k - 1] - P[i_curr][k]) 
            - beta * Q[i_curr][k] * abs(Q[i_curr][k]));
    }

    // One-spatial-step calculation of pressure
    double calculateP(int t_step, int s_step)
    {
        int i_next = (t_step + 1) % 2;
        int i_curr = t_step % 2;
        int k = s_step;

        return P[i_next][k] = P[i_curr][k] + h * (gamma * (Q[i_curr][k] - Q[i_curr][k + 1]));
    }

    void setInitP(double P_init) 
    { 
        P[0][steps - 1] = P_init;
        P[1][steps - 1] = P_init; 
    }

    std::vector<double> getLastQ() const { return std::vector<double>( {Q[1][1], Q[1][steps / 2], Q[1][steps - 1]}); };
    std::vector<double> getLastP() const { return std::vector<double>( {P[1][1], P[1][steps / 2], P[1][steps - 2]}); };

    int getSteps() const { return steps; };
    int getLength() const { return length; };
    int getId() const { return id; };
    double getGamma() const { return gamma; };
};

// Property map for vertex
class VertexProperty
{
private:
    // Also a circular buffer that describes pressure for 2 timesteps 
    double P[2];
    double h;

    // Registry object that contains incoming and outocming air flows for further adaptation
    std::vector<double> regi;

public:
    VertexProperty() : h(0), regi(0) 
    {
        P[0] = 0.0;
        P[1] = 0.0;
    };

    VertexProperty(double hh, double p) : h(hh)
    {
        P[0] = p;
        P[1] = p;
    };

    /**
     *  Adaptation of pressure inside vertex consitutes change in pressure
     *  with respect to a neighboring vertices. By Kirchhoff's 1st law, the sum of 
     *  current (pressure) entering a junction is equal to a sum of current (pressure)
     *  leaving the junction.
    */
    double adapt(int t_step, double gamma)
    {
        // Pressure calculation 
        int i_next = (t_step + 1) % 2;
        int i_curr = t_step % 2;
        double resQ;

        for(int i = 0; i < regi.size(); i++)
            resQ += regi[i];
        
        regi.clear();
        return P[i_next] = P[i_curr] + h * (gamma * resQ);
    };

    void addQ(double d) { regi.push_back(d); }
    auto getP(int t_step) const { return P[t_step % 2]; };
};

class EdgeDebug
{
private:
    EdgeProperty *edge;

    std::ofstream Q_plots_s;
    std::ofstream Q_plots_m;
    std::ofstream Q_plots_e;

    std::ofstream P_plots_s;
    std::ofstream P_plots_m;
    std::ofstream P_plots_e;

public:
    EdgeDebug(std::string name, EdgeProperty *e);
    void serialize(double time);
    int getId() const { return edge->getId(); };
};

class VertexDebug
{
private:
    size_t id;
    std::ofstream V_plots_s;

public:
    VertexDebug(std::string name) {};
};