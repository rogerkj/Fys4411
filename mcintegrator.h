#ifndef MCINTEGRATOR_H
#define MCINTEGRATOR_H

#include <armadillo>

using namespace arma;

class MCIntegrator
{
public:
    MCIntegrator();
    ~MCIntegrator();

    double runMCIntegration(double _alpha, double _beta);

private:

    void quantumforce (const mat &r , mat &qforce ,double wf);

    double waveFunction(const mat &r);
    double jastrowFactor(const mat &r);
    double localEnergy(const mat &r);
    double localEnergyAnalytical(const mat &r);

    mat r_old,r_new;
    mat qforce_old,qforce_new;

    int nDimensions;
    int charge;
    double stepLength;
    int nParticles;

    double h;
    double h2;

    long idum;

    double alpha;
    double beta;

    int nCycles;
    
    double timestep;
    double D;

    int  termiLim;

};

#endif // MCINTEGRATOR_H
