// @author Wushi Dong

// Graphene.cpp

#include "Graphene.h"

Graphene::Graphene(const Parameters &rParameters)
{
  mName = "Graphene";

  mLatticeConstant = 2.46;

  mNumWannier = 2;

  mHopping = -2.8;

  mFermiLevel = rParameters.mFermiLevelGr;

  mNumUnitCellX = 2;
  mNumUnitCellY = 4;

  mUnitCellHamiltonian.set_size(mNumWannier * mNumUnitCellX, mNumWannier
      * mNumUnitCellX); 
  SetUnitCellHamiltonian();

  mUnitCellInteractionX.set_size(mNumWannier * mNumUnitCellX, mNumWannier
      * mNumUnitCellX); 
  SetUnitCellInteractionX();

  mUnitCellInteractionY.set_size(mNumWannier * mNumUnitCellX, mNumWannier
      * mNumUnitCellX); 
  SetUnitCellInteractionY();
}


Graphene::~Graphene()
{
//  cout << "Graphene deleted" << endl;
}


//void Graphene::SetFermiLevel(const Parameters *parameters)
//{
//  mFermiLevel = parameters->mGrapheneFermiLevel;
//}


// Uses a unit cell of 4 carbon atoms
void Graphene::SetUnitCellHamiltonian()
{
  mUnitCellHamiltonian << -mFermiLevel << mHopping << 0.0 << 0.0 << endr
                       << mHopping << -mFermiLevel << mHopping << 0.0 << endr
                       << 0.0 << mHopping << -mFermiLevel << mHopping << endr
                       << 0.0 << 0.0 << mHopping << -mFermiLevel << endr;
}


void Graphene::SetUnitCellInteractionX()
{
  mUnitCellInteractionX << 0.0  << 0.0  << 0.0  << 0.0  << endr
                        << 0.0  << 0.0  << 0.0  << 0.0  << endr
                        << 0.0  << 0.0  << 0.0  << 0.0  << endr
                        << mHopping << 0.0  << 0.0  << 0.0  << endr;
}


void Graphene::SetUnitCellInteractionY()
{
  mUnitCellInteractionY << 0.0  << 0.0  << 0.0  << 0.0  << endr
                        << mHopping << 0.0  << 0.0  << 0.0  << endr
                        << 0.0  << 0.0  << 0.0  << mHopping << endr
                        << 0.0  << 0.0  << 0.0  << 0.0  << endr;
}


cx_mat Graphene::GetPrincipalLayerHamiltonian(const double phase) const
{
  cx_mat PrincipalLayerHamiltonian(mNumWannier * mNumUnitCellX * mNumUnitCellY,
      mNumWannier * mNumUnitCellX * mNumUnitCellY);

  // Used by kron() to construct the Hamiltonian matrix of the supercell
  mat f_Gr(mNumUnitCellY, mNumUnitCellY); 

  // Sets the off-diagonal entries 
  f_Gr.zeros(); 
  f_Gr.diag(1) += 1; 
  PrincipalLayerHamiltonian = kron(f_Gr, mUnitCellInteractionY); 
  PrincipalLayerHamiltonian = PrincipalLayerHamiltonian
    + PrincipalLayerHamiltonian.t();

  // Sets the diagonal entries
  f_Gr.zeros(); 
  f_Gr.diag(0) += 1; 
  PrincipalLayerHamiltonian = PrincipalLayerHamiltonian + kron(f_Gr,
      mUnitCellHamiltonian);

  // Sets the interaction with the neighbouring Cell along the periodic y axis
  PrincipalLayerHamiltonian(mNumWannier * mNumUnitCellX * mNumUnitCellY - 3, 0) 
    += mHopping * exp(constants::I * phase);
  PrincipalLayerHamiltonian(mNumWannier * mNumUnitCellX * mNumUnitCellY - 2, 3) 
    += mHopping * exp(constants::I * phase);
  PrincipalLayerHamiltonian(0, mNumWannier * mNumUnitCellX * mNumUnitCellY - 3) 
    += mHopping * exp(-constants::I * phase); 
  PrincipalLayerHamiltonian(3, mNumWannier * mNumUnitCellX * mNumUnitCellY - 2) 
    += mHopping * exp(-constants::I * phase);

  return PrincipalLayerHamiltonian;
}


// This term curently does not use the phase argument because we only include
// nearest-neighbor hoppings for graphene.
cx_mat Graphene::GetPrincipalLayerInteractionX(const double phase) const
{
  cx_mat PrincipalLayerInteractionX(mNumWannier * mNumUnitCellX * mNumUnitCellY,
      mNumWannier * mNumUnitCellX * mNumUnitCellY);

  // Used by kron() to construct the Hamiltonian matrix of the supercell
  mat f_Gr(mNumUnitCellY, mNumUnitCellY);

  f_Gr.zeros(); 
  f_Gr.diag(0) += 1;
  PrincipalLayerInteractionX = kron(f_Gr, mUnitCellInteractionX);

  return PrincipalLayerInteractionX; 
}



//// Testing
//int main()
//{
//  Graphene* pGraphene = new Graphene();
//  pGraphene->Init();
//
//  cout << "Graphene unit cell Hamiltonian (eV):" << endl;
//  cout << pGraphene->GetUnitCellHamiltonian();
//  cout << endl;
//
//  cout << "Graphene unit cell interaction along x axis (eV):" << endl;
//  cout << pGraphene->GetUnitCellInteractionX();
//  cout << endl;
//
//  cout << "Graphene unit cell interaction along y axis (eV):" << endl;
//  cout << pGraphene->GetUnitCellInteractionY();
//  cout << endl;
//
//  cout << pGraphene->GetPrincipalLayerHamiltonian(0.);
//  cout<< endl;
//  cout << pGraphene->GetPrincipalLayerInteractionX(0.);
//  cout<< endl;
//  
//  return 0;
//}
