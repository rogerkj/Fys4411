#ifndef MCINTEGRATOR_H
#define MCINTEGRATOR_H

#include <armadillo>

using namespace arma;

class MCIntegrator
{
public:
    MCIntegrator();

    void runMCIntegration();

private:
    double waveFunction(const mat &r);
    double jastrowFactor(const mat &r);
    double localEnergy(const mat &r);
    double localEnergyAnalytical(const mat &r);

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

    mat rOld;
    mat rNew;
};

#endif // MCINTEGRATOR_H
