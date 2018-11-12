// @author Wushi Dong

// FixedChargeProfile.cpp

#include "FixedChargeProfile.h"

FixedChargeProfile::~FixedChargeProfile(){}

void FixedChargeProfile::GetFixedCharge()
{
  double * charge_fixed_Gr = (double *)calloc(6, sizeof(double));
  double * charge_fixed_MoS2 = (double *)calloc(6, sizeof(double));
  double temp;
  std::ifstream infile;

  // Reads fixed charge from file
  // Graphene
  int n_scat_Gr = 10;
  std::string filename1 = mFixedChargeFileGr;
  infile.open(filename1.c_str());
  if(!infile && mRank == 0)
  {
    std::cerr<< endl << "Error when trying to open Gr fixed charge file!" << std::endl;
  }
  for(int i = 0; i < 6; ++i)
  {
    infile >> temp;
    infile >> temp;
  } 
 for(int i = 6; i < (n_scat_Gr - 1) * 6; ++i)
  {
    infile >> temp;
    infile >> temp;
    charge_fixed_Gr[i % 6] += temp;
  } 
  infile.close();

  for(int i = 0; i < 6; ++i)
    charge_fixed_Gr[i] /= (double)(n_scat_Gr - 2);

  // Displays fixed charge for graphene
  if(mRank == 0)
  {
    std::cout << "Fixed charge for Gr supercell:" << endl;
    for (int i = 0; i < 6; ++i)
      std::cout << i + 1 << "\t" << std::setprecision(7) << charge_fixed_Gr[i] << endl;
    std::cout << endl;
  }

  // MoS2
  int n_scat_MoS2 = 10;
  std::string filename2 = mFixedChargeFileMoS2;
  infile.open(filename2.c_str());
  if(!infile && mRank == 0)
  {
    std::cerr<< endl << "Error when trying to open MoS2 fixed charge file!" << std::endl;
  }
  for(int i = 0; i < 6; ++i)
  {
    infile >> temp;
    infile >> temp;
  } 
 for(int i = 6; i < (n_scat_MoS2 - 1) * 6; ++i)
  {
    infile >> temp;
    infile >> temp;
    charge_fixed_MoS2[i % 6] += temp;
  } 
  infile.close();

  for(int i = 0; i < 6; ++i)
    charge_fixed_MoS2[i] /= (double)(n_scat_MoS2 - 2);

  // Displays fixed charge for MoS2
  if(mRank == 0)
  {
    std::cout << "Fixed charge for MoS2 supercell:" << endl;
    for (int i = 0; i < 6; ++i)
      std::cout << i + 1 << "\t" << std::setprecision(7) << charge_fixed_MoS2[i] << endl;
    std::cout << endl;
  }

  // Sets fixed charge for the device
  // Graphene
  for(int i = 0; i < mNumPrincipalLayerGraphene * 6 - 1; ++i)
    mpCharge[i] = charge_fixed_Gr[i % 6];
  // S2
  mpCharge[mNumPrincipalLayerGraphene * 6] = charge_fixed_MoS2[4];
  // Gr
  for(int i = 0; i < mNumPrincipalLayerMoS2 * 6 - 1; ++i)
    mpCharge[i + mNumPrincipalLayerGraphene * 6 + 1] = charge_fixed_MoS2[i % 6];
}

//void FixedChargeProfile::Print() const
//{
//  std::cout << "Fixed Charge: " << endl;
//  for(int i = 0; i < mNGrid; ++i)
//    std::cout << i << '\t' << mpCharge[i] << endl;
//  std::cout << endl;
//}

