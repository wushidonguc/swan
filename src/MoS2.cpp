// @author Wushi Dong

// MoS2.cpp


#include "MoS2.h"

MoS2::~MoS2()
{
//  std::cout << "MoS2 deleted" << endl;
}


void MoS2::SetFermiLevel(Parameters &rParameters)
{
  const std::string doping = rParameters.mDoping;

  if(doping == "1e13")
    mFermiLevel = -2.1839 + 0.849121; // doping with electron density 1e13 cm^-2
  else if(doping == "5e13")
    mFermiLevel = -2.1839 + 0.972237;
  else if(doping == "6e13")
    mFermiLevel = -2.1839 + 0.990981;
  else if(doping == "7e13")
    mFermiLevel = -2.1839 + 1.0077;
  else if(doping == "8e13")
    mFermiLevel = -2.1839 + 1.02242;
  else if(doping == "1e14")
    mFermiLevel = -2.1839 + 1.04638; 
  else if(doping == "2e14")
    mFermiLevel = -2.1839 + 1.11189; 
  else if(doping == "3e14")
    mFermiLevel = -2.1839 + 1.15117; 
  else if(doping == "4e14")
    mFermiLevel = -2.1839 + 1.18351; 
  else if(doping == "5e14")
    mFermiLevel = -2.1839 + 1.21213; 
  else if(doping == "1e15")
    mFermiLevel = -2.1839 + 1.31826; 
  else if(doping == "1e13_low")
    mFermiLevel = -2.1839 + 0.879598; // at 4.2K
  else if(doping == "1e14_low")
    mFermiLevel = -2.1839 + 1.06333;
  else if(doping == "2e14_low")
    mFermiLevel = -2.1839 + 1.12083;
  else if(doping == "3e14_low")
    mFermiLevel = -2.1839 + 1.15578;
  else
  {
    if(mRank == 0)
    {
      std::cout << "ERROR: doping level not available!" << endl;
    }
  }

  if(mRank == 0)
  {
    std::cout << "Determined MoS2 Fermi level shift due to doping is: " <<
    mFermiLevel - (-2.1839) << " eV." << endl;
  }
}


void MoS2::SetHoppingData()
{
  std::ifstream infile;
  std::string line;
  double temp;
  int nrpt_MoS2;
  int n_r1, n_r2;
  int wan_n, wan_m;

//  std::cout << "[Material] Reading MoS2 hopping paramters from file " <<
//  mHoppingFile << " ..." << endl << endl;

//  std::cout << mHoppingFile << endl;

  infile.open(mHoppingFile.c_str());
  if(!infile && mRank == 0)
  {
    std::cerr<< endl << "Error when trying to open MoS2 hopping parameters \
      file!" << std::endl;
//    throw "Error when trying to open MoS2 hopping parameters file!";
  }
 
  // Begins to parse the file
  // Passes the first line of date and time
  getline(infile, line);  

  // Skips reading the second line of the number of the wannier functions
  infile >> temp;  

  // Reads the third line of nrpt
  infile >> nrpt_MoS2;

  //  Skips the line(s) describing the degeneracy of grid-points
  for(int i = 0; i < nrpt_MoS2 / (2 * mNumNearestNeighbors + 1) + 2; ++i)
  {
    getline(infile, line);
  }

  // Reads the parameters
  for(int i = 0; i < nrpt_MoS2 * mNumWannier * mNumWannier; ++i)
  {
    infile >> n_r1; 
    infile >> n_r2; 
    infile >> temp; 
    infile >> wan_n;
    infile >> wan_m;
    
    //  Saves the real part of the hopping parameter
    infile >> mpHoppingData[(n_r1 + mNumNearestNeighbors) * ((2
          * mNumNearestNeighbors + 1) * mNumWannier * mNumWannier) + (n_r2
          + mNumNearestNeighbors) * (mNumWannier * mNumWannier) + (wan_n - 1)
          * mNumWannier + (wan_m - 1)];

    // Shifts the Fermi level to 0 eV for the on-site terms
    if(n_r1 == 0 && n_r2 == 0 && wan_n == wan_m)
    {
      mpHoppingData[(n_r1 + mNumNearestNeighbors) * ((2 * mNumNearestNeighbors
            + 1) * mNumWannier * mNumWannier) + (n_r2 + mNumNearestNeighbors)
        * (mNumWannier * mNumWannier) +(wan_n - 1) * mNumWannier + (wan_m - 1)]
        -= mFermiLevel;
    }

    // Skips the imaginary part of hopping
    infile >> temp; 
  }

  infile.close();

}


cx_mat MoS2::GetPrincipalLayerHamiltonian(const double phase) const
{
  cx_mat PrincipalLayerHamiltonian(mNumWannier * mNumUnitCellX * mNumUnitCellY,
      mNumWannier * mNumUnitCellX * mNumUnitCellY);

  // Assign Hamiltonian for all pairs of unit cells within a principal layer
  for(int i_rep = 0; i_rep < mNumUnitCellY; ++i_rep)
  {
    for(int j_rep = 0; j_rep < mNumUnitCellY; ++j_rep)
      PrincipalLayerHamiltonian.submat(i_rep * mNumWannier * 2, j_rep
          * mNumWannier * 2, (i_rep + 1) * mNumWannier * 2 - 1, (j_rep + 1)
          * mNumWannier * 2 - 1) = ConstructMoS2Hamiltonian(0, j_rep - i_rep,
            phase);
  }	

  return PrincipalLayerHamiltonian;
}


cx_mat MoS2::GetPrincipalLayerInteractionX(const double phase) const
{
  cx_mat PrincipalLayerInteractionX(mNumWannier * mNumUnitCellX * mNumUnitCellY,
      mNumWannier * mNumUnitCellX * mNumUnitCellY);

  // Assign Hamiltonian for all pairs of unit cells in two adjacent principal 
  // layers by calling ConstructMoS2Hamiltonian(1, ...)
  for(int i_rep = 0; i_rep < mNumUnitCellY; ++i_rep)
  {
    for(int j_rep = 0; j_rep < mNumUnitCellY; ++j_rep)
      PrincipalLayerInteractionX.submat(i_rep * mNumWannier * 2, j_rep
          * mNumWannier * 2, (i_rep + 1) * mNumWannier * 2 - 1, (j_rep + 1)
          * mNumWannier * 2 - 1) = ConstructMoS2Hamiltonian(1, j_rep - i_rep,
            phase);
  }

  return PrincipalLayerInteractionX; 
}


// TODO(dongws@uchicago.edu): Remove data and directly access hopping parameters
// through mpHoppingData[]
cx_mat MoS2::ConstructMoS2Hamiltonian(const int x, const int y, const double
    phase) const
{
  // Array for storing MoS2 hopping parameters
  double data[(2 * mNumNearestNeighbors + 1)][(2 * mNumNearestNeighbors
      + 1)][mNumWannier][mNumWannier];

  // Gets hopping parameters from mpHoppingData[]
  for(int i = 0; i < 2 * mNumNearestNeighbors + 1; ++i)
    for(int j = 0; j < 2 * mNumNearestNeighbors + 1; ++j)
      for(int l = 0; l < mNumWannier; ++l)
        for(int m = 0; m < mNumWannier; ++m)
          data[i][j][l][m] = mpHoppingData[i * (2 * mNumNearestNeighbors + 1)
            * mNumWannier * mNumWannier + j * mNumWannier * mNumWannier
            + l * mNumWannier + m];

  cx_mat hamiltonian(mNumWannier * mNumUnitCellX, mNumWannier * mNumUnitCellX);

  // Assign Hamiltonian between all atom pairs
  // Mo1
  for(int i_wan = 0; i_wan < mNumWannierMo; ++i_wan)
  {
    // - Mo1
    for(int j_wan = 0; j_wan < mNumWannierMo; ++j_wan)
      hamiltonian(i_wan, j_wan) = data[0 + mNumNearestNeighbors + (x - y)][0
        + mNumNearestNeighbors + (2 * x)][i_wan][j_wan] + data[0
        + mNumNearestNeighbors + (x - (y + mNumUnitCellY))][0
        + mNumNearestNeighbors + (2 * x)][i_wan][j_wan] * exp(constants::I * phase)
        + data[0 + mNumNearestNeighbors + (x - (y - mNumUnitCellY))][0
        + mNumNearestNeighbors + (2 * x)][i_wan][j_wan] * exp(-constants::I * phase);

    // - S1
    for(int j_wan = mNumWannierMo; j_wan < mNumWannier; ++j_wan)
      hamiltonian(i_wan, j_wan) = data[0 + mNumNearestNeighbors + (x - y)][1
        + mNumNearestNeighbors + (2 * x)][i_wan][j_wan] + data[0
        + mNumNearestNeighbors + (x - (y + mNumUnitCellY))][1
        + mNumNearestNeighbors + (2 * x)][i_wan][j_wan] * exp(constants::I * phase)
        + data[0 + mNumNearestNeighbors + (x - (y - mNumUnitCellY))][1
        + mNumNearestNeighbors + (2 * x)][i_wan][j_wan] * exp(-constants::I * phase);

    // - Mo2
    for(int j_wan = 0 + mNumWannier; j_wan < mNumWannierMo + mNumWannier;
        ++j_wan)
      hamiltonian(i_wan, j_wan) = data[0 + mNumNearestNeighbors + (x - y)][1
        + mNumNearestNeighbors + (2 * x)][i_wan][j_wan - mNumWannier] + data[0
        + mNumNearestNeighbors + (x - (y + mNumUnitCellY))][1
        + mNumNearestNeighbors + (2 * x)][i_wan][j_wan - mNumWannier] * exp(constants::I
            * phase) + data[0 + mNumNearestNeighbors + (x - (y
                - mNumUnitCellY))][1 + mNumNearestNeighbors + (2
                * x)][i_wan][j_wan - mNumWannier] * exp(-constants::I * phase);

    // - S2
    for(int j_wan = mNumWannierMo + mNumWannier; j_wan < mNumWannier
        + mNumWannier; ++j_wan)
      hamiltonian(i_wan, j_wan) = data[1 + mNumNearestNeighbors + (x - y)][2
        + mNumNearestNeighbors + (2 * x)][i_wan][j_wan - mNumWannier] + data[1
        + mNumNearestNeighbors + (x - (y + mNumUnitCellY))][2
        + mNumNearestNeighbors + (2 * x)][i_wan][j_wan - mNumWannier] * exp(constants::I
            * phase) + data[1 + mNumNearestNeighbors + (x - (y
                - mNumUnitCellY))][2 + mNumNearestNeighbors + (2
                * x)][i_wan][j_wan - mNumWannier] * exp(-constants::I * phase);
  }

  // S1
  for(int i_wan = mNumWannierMo; i_wan < mNumWannier; ++i_wan)
  {
    // - Mo1
    for(int j_wan = 0; j_wan < mNumWannierMo; ++j_wan)
      hamiltonian(i_wan, j_wan) = data[0 + mNumNearestNeighbors + (x - y)][-1
        + mNumNearestNeighbors + (2 * x)][i_wan][j_wan] + data[0
        + mNumNearestNeighbors + (x - (y + mNumUnitCellY))][-1
        + mNumNearestNeighbors + (2 * x)][i_wan][j_wan] * exp(constants::I * phase)
        + data[0 + mNumNearestNeighbors + (x - (y - mNumUnitCellY))][-1
        + mNumNearestNeighbors + (2 * x)][i_wan][j_wan] * exp(-constants::I * phase);

    // - S1
    for(int j_wan = mNumWannierMo; j_wan < mNumWannier; ++j_wan)
      hamiltonian(i_wan, j_wan) = data[0 + mNumNearestNeighbors + (x - y)][0
        + mNumNearestNeighbors + (2 * x)][i_wan][j_wan] + data[0
        + mNumNearestNeighbors + (x - (y + mNumUnitCellY))][0
        + mNumNearestNeighbors + (2 * x)][i_wan][j_wan] * exp(constants::I * phase)
        + data[0 + mNumNearestNeighbors + (x - (y - mNumUnitCellY))][0
        + mNumNearestNeighbors + (2 * x)][i_wan][j_wan] * exp(-constants::I * phase);

    // - Mo2
    for(int j_wan = 0 + mNumWannier; j_wan < mNumWannierMo + mNumWannier;
        ++j_wan)
      hamiltonian(i_wan, j_wan) = data[0 + mNumNearestNeighbors + (x - y)][0
        + mNumNearestNeighbors + (2 * x)][i_wan][j_wan - mNumWannier] + data[0
        + mNumNearestNeighbors + (x - (y + mNumUnitCellY))][0
        + mNumNearestNeighbors + (2 * x)][i_wan][j_wan - mNumWannier] * exp(constants::I
            * phase) + data[0 + mNumNearestNeighbors + (x - (y
                - mNumUnitCellY))][0 + mNumNearestNeighbors + (2
                * x)][i_wan][j_wan - mNumWannier] * exp(-constants::I * phase);

    // - S2
    for(int j_wan = mNumWannierMo + mNumWannier; j_wan < mNumWannier
        + mNumWannier; ++j_wan)
      hamiltonian(i_wan, j_wan) = data[1 + mNumNearestNeighbors + (x - y)][1
        + mNumNearestNeighbors + (2 * x)][i_wan][j_wan - mNumWannier] + data[1
        + mNumNearestNeighbors + (x - (y + mNumUnitCellY))][1
        + mNumNearestNeighbors + (2 * x)][i_wan][j_wan - mNumWannier] * exp(constants::I
            * phase) + data[1 + mNumNearestNeighbors + (x - (y
                - mNumUnitCellY))][1 + mNumNearestNeighbors + (2
                * x)][i_wan][j_wan - mNumWannier] * exp(-constants::I * phase);
  }

  // Mo2
  for(int i_wan = 0 + mNumWannier; i_wan < mNumWannierMo + mNumWannier; ++i_wan)
  {
    // - Mo1
    for(int j_wan = 0; j_wan < mNumWannierMo; ++j_wan)
      hamiltonian(i_wan, j_wan) = data[0 + mNumNearestNeighbors + (x - y)][-1
        + mNumNearestNeighbors + (2 * x)][i_wan - mNumWannier][j_wan] + data[0
        + mNumNearestNeighbors + (x - (y + mNumUnitCellY))][-1
        + mNumNearestNeighbors + (2 * x)][i_wan - mNumWannier][j_wan] * exp(constants::I
            * phase) + data[0 + mNumNearestNeighbors + (x - (y
                - mNumUnitCellY))][-1 + mNumNearestNeighbors + (2 * x)][i_wan
            - mNumWannier][j_wan] * exp(-constants::I * phase);

    // - S1
    for(int j_wan = mNumWannierMo; j_wan < mNumWannier; ++j_wan)
      hamiltonian(i_wan, j_wan) = data[0 + mNumNearestNeighbors + (x - y)][0
        + mNumNearestNeighbors + (2 * x)][i_wan - mNumWannier][j_wan] + data[0
        + mNumNearestNeighbors + (x - (y + mNumUnitCellY))][0
        + mNumNearestNeighbors + (2 * x)][i_wan - mNumWannier][j_wan] * exp(constants::I
            * phase) + data[0 + mNumNearestNeighbors + (x - (y
                - mNumUnitCellY))][0 + mNumNearestNeighbors + (2 * x)][i_wan
            - mNumWannier][j_wan] * exp(-constants::I * phase);

    // - Mo2
    for(int j_wan = 0 + mNumWannier; j_wan < mNumWannierMo + mNumWannier;
        ++j_wan)
      hamiltonian(i_wan, j_wan) = data[0 + mNumNearestNeighbors + (x - y)][0
        + mNumNearestNeighbors + (2 * x)][i_wan - mNumWannier][j_wan
        - mNumWannier] + data[0 + mNumNearestNeighbors + (x - (y
              + mNumUnitCellY))][0 + mNumNearestNeighbors + (2 * x)][i_wan
        - mNumWannier][j_wan - mNumWannier] * exp(constants::I * phase) + data[0
        + mNumNearestNeighbors + (x - (y - mNumUnitCellY))][0
        + mNumNearestNeighbors + (2 * x)][i_wan - mNumWannier][j_wan
        - mNumWannier] * exp(-constants::I * phase);

    // - S2
    for(int j_wan = mNumWannierMo + mNumWannier; j_wan < mNumWannier
        + mNumWannier; ++j_wan)
      hamiltonian(i_wan, j_wan) = data[1 + mNumNearestNeighbors + (x - y)][1
        + mNumNearestNeighbors + (2 * x)][i_wan - mNumWannier][j_wan
        - mNumWannier] + data[1 + mNumNearestNeighbors + (x - (y
              + mNumUnitCellY))][1 + mNumNearestNeighbors + (2 * x)][i_wan
        - mNumWannier][j_wan - mNumWannier] * exp(constants::I * phase) + data[1
        + mNumNearestNeighbors + (x - (y - mNumUnitCellY))][1
        + mNumNearestNeighbors + (2 * x)][i_wan - mNumWannier][j_wan
        - mNumWannier] * exp(-constants::I * phase);
  }

  // S2
  for(int i_wan = mNumWannierMo + mNumWannier; i_wan < mNumWannier
      + mNumWannier; ++i_wan)
  {
    // - Mo1
    for(int j_wan = 0; j_wan < mNumWannierMo; ++j_wan)
      hamiltonian(i_wan, j_wan) = data[-1 + mNumNearestNeighbors + (x - y)][-2
        + mNumNearestNeighbors + (2 * x)][i_wan - mNumWannier][j_wan] + data[-1
        + mNumNearestNeighbors + (x - (y + mNumUnitCellY))][-2
        + mNumNearestNeighbors + (2 * x)][i_wan - mNumWannier][j_wan] * exp(constants::I
            * phase) + data[-1 + mNumNearestNeighbors + (x - (y
                - mNumUnitCellY))][-2 + mNumNearestNeighbors + (2 * x)][i_wan
            - mNumWannier][j_wan] * exp(-constants::I * phase);

    // - S1
    for(int j_wan = mNumWannierMo; j_wan < mNumWannier; ++j_wan)
      hamiltonian(i_wan, j_wan) = data[-1 + mNumNearestNeighbors + (x - y)][-1
        + mNumNearestNeighbors + (2 * x)][i_wan - mNumWannier][j_wan] + data[-1
        + mNumNearestNeighbors + (x - (y + mNumUnitCellY))][-1
        + mNumNearestNeighbors + (2 * x)][i_wan - mNumWannier][j_wan] * exp(constants::I
            * phase) + data[-1 + mNumNearestNeighbors + (x - (y
                - mNumUnitCellY))][-1 + mNumNearestNeighbors + (2 * x)][i_wan
            - mNumWannier][j_wan] * exp(-constants::I * phase);

    // - Mo2
    for(int j_wan = 0 + mNumWannier; j_wan < mNumWannierMo + mNumWannier;
        ++j_wan)
      hamiltonian(i_wan, j_wan) = data[-1 + mNumNearestNeighbors + (x - y)][-1
        + mNumNearestNeighbors + (2 * x)][i_wan - mNumWannier][j_wan
        - mNumWannier] + data[-1 + mNumNearestNeighbors + (x - (y
              + mNumUnitCellY))][-1 + mNumNearestNeighbors + (2 * x)][i_wan
        - mNumWannier][j_wan - mNumWannier] * exp(constants::I * phase) + data[-1
        + mNumNearestNeighbors + (x - (y - mNumUnitCellY))][-1
        + mNumNearestNeighbors + (2 * x)][i_wan - mNumWannier][j_wan
        - mNumWannier] * exp(-constants::I * phase);

    // - S2
    for(int j_wan = mNumWannierMo + mNumWannier; j_wan < mNumWannier
        + mNumWannier; ++j_wan)
      hamiltonian(i_wan, j_wan) = data[0 + mNumNearestNeighbors + (x - y)][0
        + mNumNearestNeighbors + (2 * x)][i_wan - mNumWannier][j_wan
        - mNumWannier] + data[0 + mNumNearestNeighbors + (x - (y
              + mNumUnitCellY))][0 + mNumNearestNeighbors + (2 * x)][i_wan
        - mNumWannier][j_wan - mNumWannier]* exp(constants::I * phase) + data[0
        + mNumNearestNeighbors + (x - (y - mNumUnitCellY))][0
        + mNumNearestNeighbors + (2 * x)][i_wan - mNumWannier][j_wan
        - mNumWannier]* exp(-constants::I * phase);
  }
  
  return hamiltonian;
}





//int main()
//{
//  Parameters parameters;
//  parameters.ParseInputFile();
//  parameters.Print();
//
//  std::cout << endl << "Setting up materials ..." << endl;
//
//  // Graphene
//  const Graphene *pGraphene = new Graphene(parameters);
//  std::cout << pGraphene->GetPrefix() << " ..." << endl;
//
//  // MoS2
//  const MoS2 *pMoS2 = new MoS2(parameters, pGraphene);
//  std::cout << pMoS2->GetPrefix() << " ..." << endl;
//  std::cout << pMoS2->GetPrincipalLayerHamiltonian(0.0) << endl;
//
//  delete pGraphene;
//  delete pMoS2;
//  
//  return 0;
//}
