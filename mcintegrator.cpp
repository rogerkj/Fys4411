#include "mcintegrator.h"
#include "lib.h"

#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

MCIntegrator:: MCIntegrator() :
    nDimensions(3),
    charge(2),
    stepLength(1.0),
    nParticles(2),
    h(0.001),
    h2(1000000),
    idum(-1),
    alpha(1.6),
    beta(1.0),
    nCycles(1000000)
{
}

void MCIntegrator::runMCIntegration() {

  double waveFunctionOld = 0;
  double waveFunctionNew = 0;

  double energySum = 0;
  double energySquaredSum = 0;

  double deltaE;

  rOld = zeros<mat>( nParticles , nDimensions );
  rNew = zeros<mat>( nParticles , nDimensions );
 
  for ( int i = 0; i < nParticles ; i ++) {
    for ( int j = 0; j < nDimensions ; j ++) {
      rOld (i , j) = stepLength * ( ran2 (&idum) - 0.5) ;
    }
  }

  rNew = rOld ; 

  for(int cycle = 0; cycle < nCycles; cycle++) {

    // Store the current value of the wave function
    waveFunctionOld = waveFunction(rOld);

    // New position to test
    for(int i = 0; i < nParticles; i++) {
      for(int j = 0; j < nDimensions; j++) {
	rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
      }

      // Recalculate the value of the wave function
      waveFunctionNew = waveFunction(rNew);

      // Check for step acceptance (if yes, update position, if no, reset position)
      if(ran2(&idum) <= (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
	for(int j = 0; j < nDimensions; j++) {
	  rOld(i,j) = rNew(i,j);
	  waveFunctionOld = waveFunctionNew;
	}
      } else {
	for(int j = 0; j < nDimensions; j++) {
	  rNew(i,j) = rOld(i,j);
	}
      }
      // update energies
      deltaE = localEnergyAnalytical(rNew);
      energySum += deltaE;
      energySquaredSum += deltaE*deltaE;
    }
  }
  double energy = energySum/(nCycles * nParticles);
  double energySquared = energySquaredSum/(nCycles * nParticles);
  cout << "Energy: " << energy << " Energy (squared sum): " << energySquared << endl;
          
}

double MCIntegrator::waveFunction(const mat &r) {

    double argument = 0;
    for(int i = 0; i < nParticles; i++) {
        double rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        argument += sqrt(rSingleParticle);
    }
    return exp(-argument * alpha) * jastrowFactor(r);
}

double MCIntegrator::jastrowFactor(const mat &r) {

    rowvec r12 = r.row(1) - r.row(0);
    double r12norm = norm(r12, 2);
    return exp(r12norm / (2 * (1 + beta * r12norm)));
}


double MCIntegrator::localEnergyAnalytical(const mat &r) {
  /*
  double lngd1 = 0;
  double lngd2 = 0;
  double lngd12 = 0;
  for(int j = 0; j < nDimensions; j++) {
    lngd1 += r(0,j) * r(0,j);
    lngd2 += r(1,j) * r(1,j);
    lngd12 = (r(1,j) - r(0,j))*(r(1,j) - r(0,j));
  }
  */
  double r1 = norm(r.row(0),2);//sqrt(lngd1);
  double r2 = norm(r.row(1),2);//sqrt(lngd2);

  double r12 = norm(r.row(1) - r.row(0),2);

  double dott = dot(r.row(0), r.row(1));  

  double div = 1 + beta*r12;
  double div2 = div*div;

  double EL1 = (alpha - charge)*(1/r1 + 1/r2) + 1/r12 - alpha*alpha;

  double EL = (alpha*(r1 + r2)/r12*(1-dott/(r1*r2)) - 1/(2*div2) - 2/r12 + 2*beta/div)/(2*div2);

  return EL1 + EL;

}

double MCIntegrator::localEnergy(const mat &r) {

    mat rPlus = zeros<mat>(nParticles, nDimensions);
    mat rMinus = zeros<mat>(nParticles, nDimensions);

    rPlus = rMinus = r;

    double waveFunctionMinus = 0;
    double waveFunctionPlus = 0;

    double waveFunctionCurrent = waveFunction(r);

    // Kinetic energy

    double kineticEnergy = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionMinus = waveFunction(rMinus);
            waveFunctionPlus = waveFunction(rPlus);
            kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
            rPlus(i,j) = r(i,j);
            rMinus(i,j) = r(i,j);
        }
    }
    kineticEnergy = 0.5 * h2 * kineticEnergy / waveFunctionCurrent;

    // Potential energy
    double potentialEnergy = 0;
    double rSingleParticle = 0;
    for(int i = 0; i < nParticles; i++) {
        rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j)*r(i,j);
        }
        potentialEnergy -= charge / sqrt(rSingleParticle);
    }
    // Contribution from electron-electron potential
    double r12 = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {
            r12 = 0;
            for(int k = 0; k < nDimensions; k++) {
                r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
            }
            potentialEnergy += 1 / sqrt(r12);
        }
    }

    return kineticEnergy + potentialEnergy;
}
