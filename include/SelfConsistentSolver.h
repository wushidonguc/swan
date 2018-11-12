// @author Wushi Dong

// SelfConsistentSolver.h
// Purpose: Self-consistent solver of the electron transport. We calculate the
// electron charge and the electrostatic potential together in a self-consistent
// fashion. The converged electrostatic potential can then be used to calculate
// the transimission properties of the device.

#ifndef SELFCONSISTENTSOLVER_H
#define SELFCONSISTENTSOLVER_H

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

class SelfConsistentSolver
{
  public:
    // Constructor
    SelfConsistentSolver(const int rank, const int nprocs, const Parameters
        &rParameters, const DeviceEdgeContact &rDeviceEdgeContact):
      mRank(rank),
      mIPhaseStart((int)((double)rParameters.mNumKPoint / (double)nprocs
            * (double)rank)),
      mIPhaseEnd((int)((double)rParameters.mNumKPoint / (double)nprocs
            * (double)(rank + 1))),
      mPhaseStep(1.0 * M_PI / (double(rParameters.mNumKPoint))),
      mChargeProfile(rParameters),
      mFixedChargeProfile(rank, rParameters),
      mPotentialProfile(rParameters),
      mChargeSolver(rParameters, rDeviceEdgeContact),
      mPoissonSolver(rank, rParameters, rDeviceEdgeContact),
      mOutput(rank),
      mNumKPoint(rParameters.mNumKPoint),
      mPotentialDamping(rParameters.mPotentialDamping),
      mChargeDamping(rParameters.mChargeDamping),
      mSelfConsistentRunMax(rParameters.mSelfConsistentRunMax),
      mPotentialConv(rParameters.mPotentialConv / constants::ROOM_TEMPERATURE
          * rParameters.mTemperature),
      mEarlyStop(rParameters.mEarlyStop){
      }

    // Destructor
    ~SelfConsistentSolver();

    // Run self consistent solver
    void Run();

  private:
    // Processor rank
    int mRank;

    // Index of the starting phase number
    int mIPhaseStart;

    // Index of the ending phase number
    int mIPhaseEnd;

    // Step between two consecutive phases
    double mPhaseStep;

    ChargeProfile mChargeProfile;

    FixedChargeProfile mFixedChargeProfile;

    PotentialProfile mPotentialProfile;

    ChargeSolver mChargeSolver;

    PoissonSolver mPoissonSolver;
    
    Output mOutput;
    
    // Simulation parameters of the self-consistent solver. For more details of their meanings, please refer to "Parameters.h".
    const int mNumKPoint;
    const double mPotentialDamping;
    const double mChargeDamping;
    const int mSelfConsistentRunMax;
    const double mPotentialConv;
    const int mEarlyStop;

};


#endif


