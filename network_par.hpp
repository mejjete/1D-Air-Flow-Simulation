#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>

enum
{
    FLOW,
    PRESSURE
};

// Property map for edge
class EdgeProperty
{
private:
    int steps;
    int id;
    int head;
    int tail;
    double alpha;
    double beta;
    double gamma;
    double h;

    // Circular buffers for flow and pressure
    std::vector<std::vector<double>> Q;
    std::vector<std::vector<double>> P;

public:
    // Edge initialization. The last parameter specifies the boundary condition for pressure
    EdgeProperty(int ID, int s, int hd, int tl, double al, double bet, double gam, double hh)
        : steps(s), id(ID), head(hd), tail(tl), alpha(al), beta(bet), gamma(gam), h(hh), 
        Q(2, std::vector<double>(s, 0.0)), P(2, std::vector<double>(s, 0.0))
    {};

    // One-spatial-step calculation of air flow
    double calculateQ(int t_step, int s_step)
    {
        int i_next = (t_step + 1) % 2;
        int i_curr = t_step % 2;
        int k = s_step;

        return Q[i_next][k] = Q[i_curr][k] + h * (alpha * (P[i_curr][k - 1] - P[i_curr][k]) 
            - beta * Q[i_curr][k] * abs(Q[i_curr][k]));
    };

    // One-spatial-step calculation of pressure
    double calculateP(int t_step, int s_step)
    {
        int i_next = (t_step + 1) % 2;
        int i_curr = t_step % 2;
        int k = s_step;

        return P[i_next][k] = P[i_curr][k] + h * (gamma * (Q[i_curr][k] - Q[i_curr][k + 1]));
    };

    // Source pressure comes from the source vertex
    void setSourceP(int t_step, double P_init) 
    {
        int i_curr = t_step % 2;
        P[i_curr][0] = P_init; 
    };

    // Target pressure comes from the target vertex
    void setTargetP(int t_step, double P_init)
    {
        int i_curr = t_step % 2;
        P[i_curr][steps - 1] = P_init;
    };

    double getLastQ(int t_step, int s_step) const 
    { 
        int i_req = t_step % 2;
        return Q[i_req][s_step]; 
    };
    
    double getLastP(int t_step, int s_step) const 
    {
        int i_req = t_step % 2;
        return P[i_req][s_step];
    };

    int getSteps() const    { return steps; };
    int getID() const       { return id; };
    int getHead() const     { return head; };
    int getTail() const     { return tail; };
    double getGamma() const { return gamma; };
};

// Property map for vertex
class VertexProperty
{
private:
    // Also a circular buffer that describes pressure for 2 timesteps 
    double P[2];
    double h;
    double gamma;
    int id;

    // Set to true if this vertex has to be updated during adaptation
    bool flags;

    /** 
     *  Each vertex has to memorise incoming and outcoming edges in order
     *  to do boundary exchange
    */
    std::vector<int> incoming;
    std::vector<int> outcoming;

public:
    VertexProperty(int ID, double hh, double p, double gm, std::vector<int> in, 
        std::vector<int> out) : h(hh), gamma(gm), id(ID), flags(false), incoming(in), outcoming(out)
    {
        P[0] = p;
        P[1] = p;
    };

    double calculateP(int t_step, double flow)
    {
        // Pressure calculation 
        int i_next = (t_step + 1) % 2;
        int i_curr = t_step % 2;

         if(flags == true)
            return P[i_next]; 
        return P[i_next] = P[i_curr] + h * (gamma * flow);
    };

    const std::vector<int> &getIncoming() const     { return incoming; };
    const std::vector<int> &getOutcoming() const    { return outcoming; };

    void setGamma(double g)     { gamma = g; };
    void setP(double Pp)        { P[0] = Pp; P[1] = Pp; };
    double getGamma() const     { return gamma; };
    auto getP(int t_step) const { return P[t_step % 2]; };
    int getID() const           { return id; };
    void setDefault()           { flags = true; };
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
    int getID() const { return edge->getID(); };
};

class VertexDebug
{
private:
    VertexProperty *vertex;
    std::ofstream V_plots;

public:
    VertexDebug(std::string name, VertexProperty *vertex);
    void serialize(double time);
    int getID() const { return vertex->getID(); };
};