
// @author Wushi Dong

// CentralRegionEdgeContact.cpp


#include "CentralRegionEdgeContact.h"


CentralRegionEdgeContact::~CentralRegionEdgeContact(){}


cx_mat CentralRegionEdgeContact::GetInteractionGrapheneS2(const double phase)
  const
{
  // Width of one supercell (same for graphene and MoS2)
  const double width = mLatticeConstantMoS2 * mNumUnitCellYMoS2;
  
  cx_mat interaction_Gr_S2(mNumWannierGraphene * mNumUnitCellXGraphene
      * mNumUnitCellYGraphene, mNumUnitCellYMoS2 * mNumWannierS2);
  interaction_Gr_S2.zeros();

  double distance;

  // position profile for Gr edge carbon atoms
  double pos_Gr[mNumWannierGraphene * mNumUnitCellXGraphene
    * mNumUnitCellYGraphene][2];

  // position profile for interfacial sulfur atoms
  double pos_s[mNumUnitCellYMoS2][2];
  
  // Sets carbon atom positions
  for(int i = 0; i < mNumWannierGraphene * mNumUnitCellXGraphene; ++i)
  {
    pos_Gr[2 * i + 3][0] = 0.0;
    pos_Gr[2 * i + 3][1] = i * mLatticeConstantGraphene;
    pos_Gr[2 * i + 2][0] = -mLatticeConstantGraphene / sqrt(3.0) * 0.5;
    pos_Gr[2 * i + 2][1] = i * mLatticeConstantGraphene + 0.5
      * mLatticeConstantGraphene;
    pos_Gr[2 * i + 1][0] = -mLatticeConstantGraphene / sqrt(3.0) * 1.5;
    pos_Gr[2 * i + 1][1] = i * mLatticeConstantGraphene + 0.5
      * mLatticeConstantGraphene;
    pos_Gr[2 * i + 0][0] = -mLatticeConstantGraphene / sqrt(3.0) * 2.0;
    pos_Gr[2 * i + 0][1] = i * mLatticeConstantGraphene;
  }
  
  // Sets sulfur atom positions
  for(int i = 0; i < mNumUnitCellYMoS2; ++i)
  {
    pos_s[i][0] = mDistanceGrS2;
    pos_s[i][1] = (double) i * mLatticeConstantMoS2 + mYShift;
  }

  // Assigns interaction parameters that exponential decays with distance
  for(int i = 0; i < mNumWannierGraphene * mNumUnitCellXGraphene
      * mNumUnitCellYGraphene; ++i)
  {
    for(int j = 0; j < mNumUnitCellYMoS2; ++j)
    {
      distance = sqrt(pow((pos_Gr[i][0] - pos_s[j][0]), 2) + pow(pos_Gr[i][1]
            - pos_s[j][1], 2));
      interaction_Gr_S2(i, j * mNumWannierS2 + mNumWannierGraphene
          * mNumUnitCellXGraphene) = mHoppingCS * exp(-distance
            / mDistanceGrS2);
    }
  }
  // With supercell above
  for(int i = 0; i < mNumWannierGraphene * mNumUnitCellXGraphene
      * mNumUnitCellYGraphene; ++i)
  {
    for(int j = 0; j < mNumUnitCellYMoS2; ++j)
    {
      distance = sqrt(pow((pos_Gr[i][0] - pos_s[j][0]), 2) + pow(pos_Gr[i][1]
            - (pos_s[j][1] + width), 2));
      interaction_Gr_S2(i, j * mNumWannierS2 + mNumWannierGraphene
          * mNumUnitCellXGraphene) += mHoppingCS * exp(-distance
            / mDistanceGrS2) * exp(constants::I * phase);
    }
  }
  // With supercell below
  for(int i = 0; i < mNumWannierGraphene * mNumUnitCellXGraphene
      * mNumUnitCellYGraphene; ++i)
  {
    for(int j = 0; j < mNumUnitCellYMoS2; ++j)
    {
      distance = sqrt(pow((pos_Gr[i][0] - pos_s[j][0]), 2) + pow(pos_Gr[i][1]
            - (pos_s[j][1] - width), 2));
      interaction_Gr_S2(i, j * mNumWannierS2 + mNumWannierGraphene
          * mNumUnitCellXGraphene) += mHoppingCS * exp(-distance
            / mDistanceGrS2) * exp(-constants::I * phase);
    }
  }

  return interaction_Gr_S2;
}


// TODO(dongws@uchicago.edu): Remove data and directly access hopping parameters
// through mpHoppingData[]
cx_mat CentralRegionEdgeContact::GetHamiltonianS2(double phase) const
{ 

  // Array for storing MoS2 hopping parameters
  double data[(2 * mNumNearestNeighborsMoS2 + 1)][(2 * mNumNearestNeighborsMoS2 + 1)][mNumWannierMoS2][mNumWannierMoS2];

  // Gets hopping parameters from mpHoppingData[]
  for(int i = 0; i < 2 * mNumNearestNeighborsMoS2 + 1; ++i)
    for(int j = 0; j < 2 * mNumNearestNeighborsMoS2 + 1; ++j)
      for(int l = 0; l < mNumWannierMoS2; ++l)
        for(int m = 0; m < mNumWannierMoS2; ++m)
          data[i][j][l][m] = mpHoppingData[i * (2 * mNumNearestNeighborsMoS2
              + 1) * mNumWannierMoS2 * mNumWannierMoS2 + j * mNumWannierMoS2
            * mNumWannierMoS2 + l * mNumWannierMoS2 + m];

  // S2 layer in the middle
  cx_mat h_s2(mNumWannierS2, mNumWannierS2);
  cx_mat int_scat_1(mNumWannierS2, mNumWannierS2);
  cx_mat int_scat_2(mNumWannierS2, mNumWannierS2);
  cx_mat PL_s2(mNumWannierS2 * mNumUnitCellYMoS2, mNumWannierS2
      * mNumUnitCellYMoS2);

  // hamiltonian for a single unit cell (phase-dependent supercell coupling
  // included)
  for(int i_wan = 0; i_wan < mNumWannierS2; ++i_wan)
  {
    for(int j_wan = 0; j_wan < mNumWannierS2; ++j_wan)
      h_s2(i_wan, j_wan) = data[0 + mNumNearestNeighborsMoS2][0
        + mNumNearestNeighborsMoS2][i_wan + mNumWannierMo][j_wan
        + mNumWannierMo] + data[0 + mNumNearestNeighborsMoS2
        - mNumUnitCellYMoS2][0 + mNumNearestNeighborsMoS2][i_wan
        + mNumWannierMo][j_wan + mNumWannierMo] * exp(constants::I * phase) + data[0
        + mNumNearestNeighborsMoS2 + mNumUnitCellYMoS2][0
        + mNumNearestNeighborsMoS2][i_wan + mNumWannierMo][j_wan
        + mNumWannierMo] * exp(-constants::I * phase);
  }

  // hamiltonian between two neighboring unit Cell in y direction
  // (phase-dependent supercell coupling included)
  for(int i_wan = 0; i_wan < mNumWannierS2; ++i_wan)
  {
    for(int j_wan = 0; j_wan < mNumWannierS2; ++j_wan)
    {
      int_scat_1(i_wan, j_wan) = data[0 + mNumNearestNeighborsMoS2 - 1][0
        + mNumNearestNeighborsMoS2][i_wan + mNumWannierMo][j_wan
        + mNumWannierMo] + data[0 + mNumNearestNeighborsMoS2
        - 1 - mNumUnitCellYMoS2][0 + mNumNearestNeighborsMoS2][i_wan
        + mNumWannierMo][j_wan + mNumWannierMo] * exp(constants::I * phase) + data[0
        + mNumNearestNeighborsMoS2 - 1 + mNumUnitCellYMoS2][0
        + mNumNearestNeighborsMoS2][i_wan + mNumWannierMo][j_wan
        + mNumWannierMo] * exp(-constants::I * phase);
      int_scat_2(i_wan, j_wan) = data[0 + mNumNearestNeighborsMoS2 - 2][0
        + mNumNearestNeighborsMoS2][i_wan + mNumWannierMo][j_wan
        + mNumWannierMo] + data[0 + mNumNearestNeighborsMoS2
        - 2 - mNumUnitCellYMoS2][0 + mNumNearestNeighborsMoS2][i_wan
        + mNumWannierMo][j_wan + mNumWannierMo] * exp(constants::I * phase) + data[0
        + mNumNearestNeighborsMoS2 - 2 + mNumUnitCellYMoS2][0
        + mNumNearestNeighborsMoS2][i_wan + mNumWannierMo][j_wan
        + mNumWannierMo] * exp(-constants::I * phase);
    }
  }

  // Hamiltonian for a supercell
  mat f_scat(mNumUnitCellYMoS2, mNumUnitCellYMoS2);
  f_scat.zeros();
  f_scat.diag(1) += 1;
  PL_s2 = kron(f_scat, int_scat_1);
  f_scat.zeros();
  f_scat.diag(2) += 1;
  PL_s2 = PL_s2 + kron(f_scat, int_scat_2);
  PL_s2 = PL_s2 + PL_s2.t();
  f_scat.zeros();
  f_scat.diag(0) += 1;
  PL_s2 = PL_s2 + kron(f_scat,h_s2);

  return PL_s2;
}


cx_mat CentralRegionEdgeContact::GetInteractionS2MoS2(const double phase) const
{

  cx_mat INT_s2_MoS2(mNumWannierS2 * mNumUnitCellYMoS2, mNumWannierMoS2
      * mNumUnitCellXMoS2 * mNumUnitCellYMoS2);
  for(int i_rep = 0; i_rep < mNumUnitCellYMoS2; ++i_rep)
  {
    for(int j_rep = 0; j_rep < mNumUnitCellYMoS2; ++j_rep)
      INT_s2_MoS2.submat(i_rep * mNumWannierS2, j_rep * mNumWannierMoS2
          * mNumUnitCellXMoS2, (i_rep + 1) * mNumWannierS2 - 1, (j_rep + 1)
          * mNumWannierMoS2 * mNumUnitCellXMoS2 - 1)
        = GetS2MoS2Interaction(j_rep - i_rep, phase);
  }

  return INT_s2_MoS2;
}


// TODO(dongws@uchicago.edu): Remove data and directly access hopping parameters
// through mpHoppingData[]
cx_mat CentralRegionEdgeContact::GetS2MoS2Interaction(const int y, const double
    phase) const
{
  cx_mat hamiltonian(mNumWannierS2, mNumWannierMoS2 * mNumUnitCellXMoS2);
  
  // Array for storing MoS2 hopping parameters
  double data[(2 * mNumNearestNeighborsMoS2 + 1)][(2 * mNumNearestNeighborsMoS2
      + 1)][mNumWannierMoS2][mNumWannierMoS2];
  
  // Gets hopping parameters from mpHoppingData[]
  for(int i = 0; i < 2 * mNumNearestNeighborsMoS2 + 1; ++i)
    for(int j = 0; j < 2 * mNumNearestNeighborsMoS2 + 1; ++j)
      for(int l = 0; l < mNumWannierMoS2; ++l)
        for(int m = 0; m < mNumWannierMoS2; ++m)
          data[i][j][l][m] = mpHoppingData[i * (2 * mNumNearestNeighborsMoS2
              + 1) * mNumWannierMoS2 * mNumWannierMoS2 + j * mNumWannierMoS2
            * mNumWannierMoS2 + l * mNumWannierMoS2 + m];

  // Assign Hamiltonian between all atom pairs
  // S
  for(int i_wan = 0; i_wan < mNumWannierS2; ++i_wan)
  {
    // - Mo1
    for(int j_wan = 0; j_wan < mNumWannierMo; ++j_wan)
      hamiltonian(i_wan, j_wan) = data[0 + mNumNearestNeighborsMoS2 + (- y)][0
        + mNumNearestNeighborsMoS2][i_wan + mNumWannierMo][j_wan] + data[0
        + mNumNearestNeighborsMoS2 + (-(y + mNumUnitCellYMoS2))][0
        + mNumNearestNeighborsMoS2][i_wan + mNumWannierMo][j_wan] * exp(constants::I
            * phase) + data[0 + mNumNearestNeighborsMoS2 + (-(y
                - mNumUnitCellYMoS2))][0 + mNumNearestNeighborsMoS2][i_wan
            + mNumWannierMo][j_wan] * exp(-constants::I * phase);
    // - S1
    for(int j_wan = mNumWannierMo; j_wan < mNumWannierMoS2; ++j_wan)
      hamiltonian(i_wan, j_wan) = data[0 + mNumNearestNeighborsMoS2 + (- y)][1
        + mNumNearestNeighborsMoS2][i_wan + mNumWannierMo][j_wan] + data[0
        + mNumNearestNeighborsMoS2 + (-(y + mNumUnitCellYMoS2))][1
        + mNumNearestNeighborsMoS2][i_wan + mNumWannierMo][j_wan] * exp(constants::I
            * phase) + data[0 + mNumNearestNeighborsMoS2 + (-(y
                - mNumUnitCellYMoS2))][1 + mNumNearestNeighborsMoS2][i_wan
            + mNumWannierMo][j_wan] * exp(-constants::I * phase);
    // - Mo2
    for(int j_wan = 0 + mNumWannierMoS2; j_wan < mNumWannierMo
        + mNumWannierMoS2; ++j_wan)
      hamiltonian(i_wan, j_wan) = data[0 + mNumNearestNeighborsMoS2 + (- y)][1
        + mNumNearestNeighborsMoS2][i_wan + mNumWannierMo][j_wan
        - mNumWannierMoS2] + data[0 + mNumNearestNeighborsMoS2 + (-(y
              + mNumUnitCellYMoS2))][1 + mNumNearestNeighborsMoS2][i_wan
        + mNumWannierMo][j_wan - mNumWannierMoS2] * exp(constants::I * phase) + data[0
        + mNumNearestNeighborsMoS2 + (-(y - mNumUnitCellYMoS2))][1
        + mNumNearestNeighborsMoS2][i_wan + mNumWannierMo][j_wan
        - mNumWannierMoS2] * exp(-constants::I * phase);
    // - S2
    for(int j_wan = mNumWannierMo + mNumWannierMoS2; j_wan < mNumWannierMoS2
        + mNumWannierMoS2; ++j_wan)
      hamiltonian(i_wan, j_wan) = data[1 + mNumNearestNeighborsMoS2 + (- y)][2
        + mNumNearestNeighborsMoS2][i_wan + mNumWannierMo][j_wan
        - mNumWannierMoS2] + data[1 + mNumNearestNeighborsMoS2 + (-(y
              + mNumUnitCellYMoS2))][2 + mNumNearestNeighborsMoS2][i_wan
        + mNumWannierMo][j_wan - mNumWannierMoS2] * exp(constants::I * phase) + data[1
        + mNumNearestNeighborsMoS2 + (-(y - mNumUnitCellYMoS2))][2
        + mNumNearestNeighborsMoS2][i_wan + mNumWannierMo][j_wan
        - mNumWannierMoS2] * exp(-constants::I * phase);
  }
  return hamiltonian;
}


//// Testing
//int main()
//{
//  // Input parameters
//  Parameters parameters;
//  parameters.ParseInputFile();
//  parameters.Print();
//
//  cout << endl << "Setting up materials ..." << endl;
//
//  // Graphene
//  const Graphene *pGraphene = new Graphene(parameters);
//  cout << pGraphene->GetPrefix() << " ..." << endl;
//
//  // MoS2
//  const MoS2 *pMoS2 = new MoS2(parameters, pGraphene);
//  cout << pMoS2->GetPrefix() << " ..." << endl;
//
//  // Set up edge contact device central region
//  CentralRegionEdgeContact central_region_edge_contact(parameters, pGraphene,
//  pMoS2);
//
////  cout << "GetInteractionGrapheneS2" << endl;
////  cout << central_region_edge_contact.GetInteractionGrapheneS2(0.) << endl;
////  cout << endl;
////  cout << "GetHamiltonianS2" << endl;
////  cout << central_region_edge_contact.GetHamiltonianS2(0.) << endl;
////  cout << endl;
//  cout << "GetInteractionS2MoS2" << endl;
//  cout << central_region_edge_contact.GetInteractionS2MoS2(0.) << endl;
//  cout << endl;
//
//  delete pGraphene;
//  delete pMoS2;
//
//  return 0;
//}
