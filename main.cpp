#include "mcintegrator.h"

#include <iostream>

using namespace std;

int main (int argc, char* argv[])
{

  MCIntegrator *integrator = new MCIntegrator();
  integrator->runMCIntegration();

  return(0);

}
