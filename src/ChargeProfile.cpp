// @author Wushi Dong

// ChargeProfile.cpp

#include "ChargeProfile.h"

ChargeProfile::~ChargeProfile(){}

void ChargeProfile::SetCharge(double *charge)
{
  for(int i = 0; i < mNumGrid; ++i)
    mpCharge[i] = charge[i];
}

void ChargeProfile::Print() const
{
  cout << "Charge: " << endl;
  for(int i = 0; i < mNumGrid; ++i)
    cout << i + 1 << '\t' << mpCharge[i] << endl;
  cout << endl;
}

