// @author Wushi Dong

// Parameters.h
// Purpose: Contains all the parameters of the simulation and is responsible for
// their I/O

#include <iostream>
#include <string>
#include <fstream>

#include "header.h"

#ifndef PARAMETERS_H
#define PARAMETERS_H

/* Parameters class */

class Parameters
{
  
  public:
    // Constructor
    Parameters(); 

    // Destructor
    ~Parameters(); 

    // Read parameters from file
    void ParseInputFile();

    // Print parameters
    void Print() const;

    // Simulation temperature (K)
    double mTemperature; 

    // Number of principal layers of graphene in the central region
    int mNumPrincipalLayerGraphene;

    // Number of principal layers of MoS2 in the central region
    int mNumPrincipalLayerMoS2;

    double mFermiLevelGr;
    
    // MoS2 doping level
    std::string mDoping; 

    // Distance between graphene and MoS2
    double mDistanceGrS2;

    // Hopping between pz orbital of graphene edge carbon and px orbital of MoS2 
    // edge sulfur
    double mHoppingCS;

    // Energy minimum, Energy maximum, Energy step (eV)
    double mEnergyMin, mEnergyMax, mEnergyStep; 

//    // Dense energy range for transmission and LDOS calculations
//    double mEnergyMinDense, mEnergyMaxDense, mEnergyStepDense; 

    // Number of k-point samplings in the irreducible BZ
    int mNumKPoint; 

    // Voltage bias (V)
    double mVoltageBias; 

    // Convergence criteria for calculating lead surface Green's functions (eV)
    double mSurfaceGreensFunctionConv; 

    // Potential convergence (V)
    double mPotentialConv; 

    // Maximum number of runs for self-consistent calculation
    int mSelfConsistentRunMax; 

    // Maximum number of runs for poisson solver
    int mPoissonRunMax; 

    // Potential damping factor in self-consistent calculation
    double mPotentialDamping; 

    // Charge damping factor in self-consistent calculation
    double mChargeDamping; 

    // Potential damping factor in poisson solver
    double mPotentialDampingPoisson; 
 
    // Stop the iteration if convergence does not improve in this number of
    // consecutive runs
    int mEarlyStop; 

    // Potential initialization method
    std::string mPotentialInitializationMethod; 

    // File containing hopping parameters
    std::string mHoppingFile;
    
    // File containing graphene fixed charge profile
    std::string mFixedChargeFileGr;

    // File containing MoS2 fixed charge profile
    std::string mFixedChargeFileMoS2;
   
};

#endif

