// @author Wushi Dong

// PotentialProfile.h
// Purpose: Electrostatic potential profile class.

#ifndef POTENTIALPROFILE_H
#define POTENTIALPROFILE_H

#include "Parameters.h"
//#include "DeviceEdgeContact.h"

class PotentialProfile
{
  public:  
    // Constructor    
    PotentialProfile(const Parameters &rParameters):
      mNumGrid(6 * (rParameters.mNumPrincipalLayerGraphene
            + rParameters.mNumPrincipalLayerMoS2)),
      mNumGridGraphene(6 * rParameters.mNumPrincipalLayerGraphene),
      mPotentialInitializationMethod(rParameters.mPotentialInitializationMethod),
      mVoltageBias(rParameters.mVoltageBias){
        mpPotential = (double *)calloc(mNumGrid, sizeof(double));
        Init();
      }
    
    // Constructor    
    PotentialProfile(const Parameters &rParameters, std::string 
        potentialInitializationMethod):
      mNumGrid(6 * (rParameters.mNumPrincipalLayerGraphene
            + rParameters.mNumPrincipalLayerMoS2)),
      mNumGridGraphene(6 * rParameters.mNumPrincipalLayerGraphene),
      mPotentialInitializationMethod(potentialInitializationMethod),
      mVoltageBias(rParameters.mVoltageBias){
        mpPotential = (double *)calloc(mNumGrid, sizeof(double));
        Init();
      }

    // Destructor
    ~PotentialProfile();
    
    void SetPotential(double Potential[]);
    double *GetPotential() const;
    void Print() const;
    const int GetNumGrid() const;

  private:
    int mNumGrid;
    int mNumGridGraphene;
    double *mpPotential;
    std::string mPotentialInitializationMethod;
    double mVoltageBias;

    void Init();
};

inline
double *PotentialProfile::GetPotential() const
{
  return mpPotential;
}

inline
const int PotentialProfile::GetNumGrid() const
{
  return mNumGrid;
}

#endif

