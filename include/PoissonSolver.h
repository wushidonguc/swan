// @author Wushi Dong

// PoissonSolver.h
// Purpose: Given certain charge profile, calculates potential profile from 
// Poisson equation using Newton-Raphson iteration method.

#ifndef POISSONSOLVER_H
#define POISSONSOLVER_H

#include <iostream>
#include <iomanip>
#include "mpi.h"

#include "Parameters.h"
#include "Material.h"
#include "Graphene.h"
#include "MoS2.h"
#include "ChargeProfile.h"
#include "FixedChargeProfile.h"
#include "PotentialProfile.h"

class PoissonSolver
{
  public:
    // Constructor    
    PoissonSolver(int rank, const Parameters &rParameters, const
        DeviceEdgeContact &rDeviceEdgeContact):
      mRank(rank),
      mNumPrincipalLayerGraphene(rDeviceEdgeContact.GetNumPrincipalLayerGraphene()),
      mLatticeConstantGraphene(rDeviceEdgeContact.GetGrapheneLatticeConstant()),
      mNumPrincipalLayerMoS2(rDeviceEdgeContact.GetNumPrincipalLayerMoS2()),
      mLatticeConstantMoS2(rDeviceEdgeContact.GetMoS2LatticeConstant()),
      mKT(constants::BOLTZMANN_CONSTANT * rParameters.mTemperature),
      mDistanceGrS2(rParameters.mDistanceGrS2),
      mPotentialDampingPoisson(rParameters.mPotentialDampingPoisson),
      mPoissonRunMax(rParameters.mPoissonRunMax),
      mPotentialConv(rParameters.mPotentialConv / constants::ROOM_TEMPERATURE
          * rParameters.mTemperature){
        DIL_Graphene = 6.9;
        DIL_MoS2 = 4.0;
      }

    // Destructor
    ~PoissonSolver();

    // Run the poisson solver
    void Run(PotentialProfile potential_profile, const ChargeProfile
        charge_profile, const FixedChargeProfile fixed_charge_profile, bool
        &rIsInnerConv);

  private:
    // Rank
    int mRank;

    int mNumPrincipalLayerGraphene;
    double mLatticeConstantGraphene;
    int mNumPrincipalLayerMoS2;
    double mLatticeConstantMoS2;

    // Boltzman constant * Temperature
    double mKT;
    // Dielectric constant of graphene
    double DIL_Graphene;
    // Dielectric constant of MoS2
    double DIL_MoS2;
    // Distance between graphene and S2
    double mDistanceGrS2;
    // Potential damping factor in poisson solver
    double mPotentialDampingPoisson;
    // Maximum number of runs for poisson solver 
    int mPoissonRunMax;
    // Potential convergence (V)
    double mPotentialConv;
};


#endif
