// @author Wushi Dong

// Utils.cpp

#include "Utils.h"

double Utils::GetFermiDistribution(const double energy, const double 
    fermi_level, const double kT) const
{
  return(1.0 / (1.0 + exp((energy - fermi_level) / kT)));
}


