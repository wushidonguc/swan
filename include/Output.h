// @author Wushi Dong

// Output.h
// Purpose: Writes simulation results to file.

#include <iomanip>
#include <fstream>
#include "mpi.h"

#include "ChargeProfile.h"
#include "FixedChargeProfile.h"
#include "PotentialProfile.h"
#include "EnergyRange.h"

#ifndef OUTPUT_H
#define OUTPUT_H

class Output
{
  public:
    // Constructor
    Output(const int rank):
      mRank(rank){}

    // Destructor
    ~Output();

    // Write simulation results to file
    void WriteChargeToFile(const ChargeProfile, const FixedChargeProfile, const
        bool &) const;
    void WritePotentialToFile(const PotentialProfile, const bool &) const;
    void WriteTransmissionToFile(double *, EnergyRange) const;
    void WriteLDOSToFile(double *, int, EnergyRange) const;

  private:
    const int mRank;
};


#endif

