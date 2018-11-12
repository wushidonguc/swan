// @author Wushi Dong

// EnergyRange.h
// Purpose: Class of energy range containing min, max, spacing, and number of 
// point

#ifndef ENERGYRANGE_H
#define ENERGYRANGE_H

#include "OneDimRange.h"
#include "Parameters.h"

/* Potential profile class */

class EnergyRange : public OneDimRange
{
  public:
    // Constructor
    EnergyRange(double start, double end, int n): OneDimRange(start, end, n){}
    EnergyRange(double start, double end, double spacing): OneDimRange(start,
        end, spacing){}    
    EnergyRange(Parameters parameters): OneDimRange(parameters.mEnergyMin,
        parameters.mEnergyMax, parameters.mEnergyStep){}
    
    // Destructor
    ~EnergyRange();
};

#endif

