#include "network.hpp"

EdgeDebug::EdgeDebug(std::string name, EdgeProperty *e) : edge(e)
{
    Q_plots_s.open(name + "_Qplots_s.txt");
    Q_plots_m.open(name + "_Qplots_m.txt");
    Q_plots_e.open(name + "_Qplots_e.txt");

    P_plots_s.open(name + "_Pplots_s.txt");
    P_plots_m.open(name + "_Pplots_m.txt");
    P_plots_e.open(name + "_Pplots_e.txt");
}

void EdgeDebug::serialize(double time)
{
    auto lastQ = edge->getLastQ();
    auto lastP = edge->getLastP();

    Q_plots_s << time << " " << lastQ[0] << "\n";
    Q_plots_m << time << " " << lastQ[1] << "\n";
    Q_plots_e << time << " " << lastQ[2] << "\n";

    P_plots_s << time << " " << lastP[0] << "\n";
    P_plots_m << time << " " << lastP[1] << "\n";
    P_plots_e << time << " " << lastP[2] << "\n";
}

VertexDebug::VertexDebug(std::string name, VertexProperty *vert) : vertex(vert)
{
    V_plots.open(name + "_Pplots.txt");
}

void VertexDebug::serialize(double time)
{
    V_plots << time << " " << vertex->getP(0) << "\n";
}