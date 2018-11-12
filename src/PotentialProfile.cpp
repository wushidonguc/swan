//

#include "PotentialProfile.h"

PotentialProfile::~PotentialProfile(){}

void PotentialProfile::Init()
{
  // Use linear drop as trial potential
  if(mPotentialInitializationMethod == "linear")
  {
    for(int i = 0; i < mNumGrid; ++i)
      mpPotential[i] = -0.5 * mVoltageBias + (mVoltageBias / (double)(mNumGrid + 1)) * (double)(i + 1);
  }

  // Mannually set trial potential due to screening of interface charge potential
  else if(mPotentialInitializationMethod == "mannual")
  {
    mpPotential[mNumGridGraphene] += -1.0;
    mpPotential[mNumGridGraphene - 1] += -0.5;
    mpPotential[mNumGridGraphene + 1] += -0.5;
    mpPotential[mNumGridGraphene - 2] += -0.05;
    mpPotential[mNumGridGraphene + 2] += -0.05;
  }
  
  // Use potential from last run as input for the current calculation
  else if(mPotentialInitializationMethod == "last")
  {
    std::ifstream infile;
    double temp;
    
    infile.open("last/potential.txt");
    if(!infile)
    {
      std::cerr<< "[ERROR] Cannot find last potential file!" << std::endl;
    }
    
    for(int i = 0; i < mNumGrid; ++i)
    {
      infile >> temp;
      infile >> mpPotential[i];
    }
    
    infile.close();
  }

  // Use converged potential as input for the current calculation
  else if(mPotentialInitializationMethod == "converged")
  {
    std::ifstream infile;
    double temp;
    
    infile.open("potential_conv.txt"); 
    if(!infile)
    {
      std::cerr<< "[ERROR] Cannot find converged potential file!" << std::endl;
    }
    
    for(int i = 0; i < mNumGrid; ++i)
    {
      infile >> temp;
      infile >> mpPotential[i];
    }
    
    infile.close();
  }

  else
  {
    {
      std::cout << "[Warning] The following potential initialization method is \
        not available:" << std::endl << mPotentialInitializationMethod <<
        std::endl << "Initialize potential as all zeros." << std::endl;
    }
  }
};


void PotentialProfile::SetPotential(double *Potential)
{
  for(int i = 0; i < mNumGrid; ++i)
    mpPotential[i] = Potential[i];
};


void PotentialProfile::Print() const
{
//  std::cout << "Potential: " << std::endl;
  for(int i = 0; i < mNumGrid; ++i)
    std::cout << i + 1 << '\t' << mpPotential[i] << std::endl;
  std::cout << std::endl;
};


//// Testing
//int main()
//{
//  Parameters parameters;
//  PotentialProfile potential_profile(parameters);
//  potential_profile.Print();
//
//  const int n_grid2 = potential_profile.GetNumGrid();
//  std::cout << n_grid2 << std::endl;
//  double potential[n_grid2];
//  for(int i = 0; i < n_grid2; ++i)
//    potential[i] = 0.1 * double(i);
//  
//  potential_profile.SetPotential(potential);
//  potential_profile.Print();
//
//  PotentialProfile potential_profile2(parameters, "converged");
//  potential_profile2.Print();
//  
//
//  return 0;
//}

