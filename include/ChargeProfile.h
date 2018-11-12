// @author Wushi Dong

// ChargeProfile.h
// Purpose: Charge profile

#ifndef CHARGEPROFILE_H
#define CHARGEPROFILE_H

#include "Parameters.h"
#include "DeviceEdgeContact.h"

class ChargeProfile
{
  public:
    // Constructor
    ChargeProfile(int n_grid):
      mNumGrid(n_grid){
        mpCharge = (double *)calloc(mNumGrid, sizeof(double));
      }
    
    ChargeProfile(const Parameters &rParameters):
      mNumGrid(6 * (rParameters.mNumPrincipalLayerGraphene
            + rParameters.mNumPrincipalLayerMoS2)){
        mpCharge = (double *)calloc(mNumGrid, sizeof(double));
      }
      
    ChargeProfile(const DeviceEdgeContact &rDeviceEdgeContact):
      mNumGrid(6 * (rDeviceEdgeContact.GetNumPrincipalLayerGraphene()
            + rDeviceEdgeContact.GetNumPrincipalLayerMoS2())){
        mpCharge = (double *)calloc(mNumGrid, sizeof(double));
      }

    // Destructor
    ~ChargeProfile();
    
    const int GetNumGrid() const;

    void SetCharge(double Charge[]);

    double *GetCharge() const;

    void Print() const;

  protected:
    // Number of grid point
    int mNumGrid;

    // Charge profile
    double *mpCharge;
};


inline
const int ChargeProfile::GetNumGrid() const
{
  return mNumGrid;
}

inline
double *ChargeProfile::GetCharge() const
{
  return mpCharge;
}

#endif

