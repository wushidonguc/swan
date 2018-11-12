// @author Wushi Dong

// TransportSolver.h
// Purpose: Given the converged electrostatic potential profile, solves for the 
// transmission spectrum using the Landauer-Büttiker formula as well as the 
// Local Density of States (LDOS) of the central region.

#ifndef TRANSPORTSOLVER_H
#define TRANSPORTSOLVER_H

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "mpi.h"

#include "header.h"
#include "Parameters.h"
#include "DeviceEdgeContact.h"
#include "ChargeProfile.h"
#include "FixedChargeProfile.h"
#include "PotentialProfile.h"
#include "ChargeSolver.h"
#include "PoissonSolver.h"
#include "Output.h"

class TransportSolver
{
  public:
    // Constructor
    TransportSolver(const int rank, const int nprocs, const Parameters
        &rParameters, const DeviceEdgeContact &rDeviceEdgeContact):
      mRank(rank),
      mNumKPoint(rParameters.mNumKPoint),
      mIPhaseStart((int)((double)rParameters.mNumKPoint / (double)nprocs
            * (double)rank)),
      mIPhaseEnd((int)((double)rParameters.mNumKPoint / (double)nprocs
            * (double)(rank + 1))),
      mPhaseStep(1.0 * M_PI / (double(rParameters.mNumKPoint))),
      mEnergyRangeDense(rParameters.mEnergyMin, rParameters.mEnergyMax, 0.1
          * rParameters.mEnergyStep),
      mChargeSolver(rParameters, rDeviceEdgeContact),
      // Initialize mConvergedPotentialProfile with converged potential
      mConvergedPotentialProfile(rParameters, "converged"),
      mOutput(rank){
      }

    // Destructor
    ~TransportSolver();

    // Run transport solver to obtain transmission and LDOS
    void Run();

  private:
    // Processor rank
    int mRank;

    // Number of k-point samplings    
    const int mNumKPoint;

    // Index of the starting phase number
    int mIPhaseStart;

    // Index of the ending phase number
    int mIPhaseEnd;

    // Step between two consecutive phases
    double mPhaseStep;

    // Dense energe range
    EnergyRange mEnergyRangeDense;

    ChargeSolver mChargeSolver;
    
    PotentialProfile mConvergedPotentialProfile;

    Output mOutput;

};


#endif


