#include "mcintegrator.h"

#include <iostream>

using namespace std;

int main (int argc, char* argv[])
{

  MCIntegrator *integrator = new MCIntegrator();
  double res = integrator->runMCIntegration(1.8,1.2);

  delete integrator;
  
  return(0);

}
