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
    alpha(2.0),
    beta(1.2),
    nCycles(1000000),
    timestep(0.01),
    D(0.5),
    termiLim(100000)
{

  r_old  = zeros<mat>( nParticles , nDimensions );
  r_new  = zeros<mat>( nParticles , nDimensions );

  qforce_old  = zeros<mat>( nParticles , nDimensions );
  qforce_new  = zeros<mat>( nParticles , nDimensions );

}

MCIntegrator::~MCIntegrator() {


}

double MCIntegrator::runMCIntegration(double _alpha, double _beta) {

  alpha = _alpha;
  beta = _beta;

  double energySum = 0;
  double energySquaredSum = 0;

  double deltaE = 0.0;
 

  //Finne random start til r
  r_old.randn();
  r_old *= sqrt(timestep);

  //Finne startverdi av bølgefunskjon og quantumforce
  double wfold = waveFunction (r_old) ;
  quantumforce(r_old , qforce_old , wfold ) ;
   
  r_new = r_old; 

  //Antall Cycles
  for(int cycle = 0; cycle < nCycles; cycle++) {

    // New position to test
    for(int i = 0; i < nParticles; i++) {
      for (int j =0; j < nDimensions ; j ++) {

	  vec v = zeros<vec>(1);
	  v.randn();

	  r_new (i , j) = r_old(i, j) + v(0) * sqrt(timestep) + qforce_old(i, j) * timestep *D;
	
      }
        

      //Oppdatere matrisa
      for ( int k = 0; k < nParticles ; k++) {
	if ( k != i ) {
	  for ( int j =0; j < nDimensions ; j ++) {
	    r_new (k, j) = r_old( k, j);
	  }
	}
      }     

      //Finne bølgefunsjon og quantumforce
      double wfnew = waveFunction (r_new) ;
      quantumforce ( r_new , qforce_new , wfnew) ;
      
      //Finne greensfunkjsonen
      double greensfunction = 0.0 ;
      for ( int j = 0; j < nDimensions ; j ++) {
	greensfunction += 0.5 *( qforce_old( i, j) + qforce_new( i, j)) *
	  (D*timestep * 0.5 * (  qforce_old( i, j) - qforce_new( i, j)) - r_new(i, j) + r_old(i, j)) ;
      }

      greensfunction = exp ( greensfunction ) ;

      //Random test
      if( ran2(&idum) <= greensfunction * wfnew * wfnew / wfold / wfold ) {
	//Ny posisjon
	for ( int j =0; j < nDimensions ; j ++) {
	  r_old( i, j) = r_new(i, j) ;
	  qforce_old(i, j) = qforce_new(i, j) ;
	}
	
	wfold = wfnew ;
      } 
      /*
      else {      
	for(int j = 0; j < nDimensions; j++) {
	  r_new(i,j) = r_old(i,j);
	}
      }
      */

      //Kutte for terminalisering
      if(cycle >=  termiLim) {
      
	// update energies
	deltaE = localEnergyAnalytical(r_new);
	energySum += deltaE;
	energySquaredSum += deltaE*deltaE;
      }

    }
  }
  //Finne og plotte resultat
  double energy = energySum/((nCycles -  termiLim) * nParticles);
  double energySquared = energySquaredSum/((nCycles -  termiLim) * nParticles);
  cout << "Energy: " << energy << " Energy (squared sum): " << energySquared << endl;
  
  return energy;
        
}

//Finner verdi for quantum force. Numerisk versjon.
void MCIntegrator::quantumforce (const mat &r , mat &qforce ,double wf ) {

  double wfminus , wfplus ;
  mat r_plus  = zeros<mat>( nParticles , nDimensions );
  mat r_minus = zeros<mat>( nParticles , nDimensions );

  r_plus = r_minus = r;

  // compute the firstderivative
  for ( int i = 0; i <  nParticles ; i ++) {
    for ( int j = 0; j <  nDimensions ; j ++) {
      r_plus(i, j) =  r(i, j) + h;
      r_minus(i, j) = r(i, j) - h;

      wfminus =  waveFunction(r_minus);
      wfplus =  waveFunction(r_plus);

      qforce(i, j) = (wfplus - wfminus) * 2.0 / wf / (2.0 * h);
    
      r_plus(i, j) = r(i, j) ;
      r_minus( i, j) = r(i, j) ;
    }
  }
} // end of quantum force function


//Bølgefunskjonen 
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

//Jastrow faktor
double MCIntegrator::jastrowFactor(const mat &r) {

    rowvec r12 = r.row(1) - r.row(0);
    double r12norm = norm(r12, 2);
    return exp(r12norm / (2 * (1 + beta * r12norm)));
}

//Local energy analytical version
double MCIntegrator::localEnergyAnalytical(const mat &r) {

  double r1 = norm(r.row(0),2);
  double r2 = norm(r.row(1),2);

  double r12 = norm(r.row(1) - r.row(0),2);

  double dott = dot(r.row(0), r.row(1));  

  double div = 1 + beta*r12;
  double div2 = div*div;

  double EL1 = (alpha - charge)*(1/r1 + 1/r2) + 1/r12 - alpha*alpha;

  double EL = (alpha*(r1 + r2)/r12*(1-dott/(r1*r2)) - 1/(2*div2) - 2/r12 + 2*beta/div)/(2*div2);

  return EL1 + EL;

}

//Local energy numerisk versjon
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
