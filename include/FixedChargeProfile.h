// @author Wushi Dong

// FixedChargeProfile.h
// Purpose: 

#ifndef FIXEDCHARGEPROFILE_H
#define FIXEDCHARGEPROFILE_H

#include <iostream>
#include <iomanip>
#include <string>

#include "Parameters.h"
#include "ChargeProfile.h"

class FixedChargeProfile : public ChargeProfile
{
  public:
    // Constructor
    FixedChargeProfile(const int rank, const Parameters &rParameters):
      ChargeProfile(rParameters),
      mRank(rank),
      mFixedChargeFileGr(rParameters.mFixedChargeFileGr),
      mNumPrincipalLayerGraphene(rParameters.mNumPrincipalLayerGraphene),
      mFixedChargeFileMoS2(rParameters.mFixedChargeFileMoS2),
      mNumPrincipalLayerMoS2(rParameters.mNumPrincipalLayerMoS2)      
      {
        std::cout << "Reading the fixed charges" << endl;
        GetFixedCharge();
      }
      
    // Destructor
    ~FixedChargeProfile();

//    virtual void Print() const;

  private:
    // Rank
    const int mRank;
    
    std::string mFixedChargeFileGr;
    int mNumPrincipalLayerGraphene;

    std::string mFixedChargeFileMoS2;
    int mNumPrincipalLayerMoS2;

    void GetFixedCharge();
};

#endif
