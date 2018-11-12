// @author Wushi Dong

// MoS2.h
// Purpose: Class that contains information for material MoS2

#ifndef MOS2_H
#define MOS2_H

#include <iostream>
#include <fstream>
//#include <map>
#include <string>
#include "mpi.h"

#include "header.h"
#include "Parameters.h"
#include "Material.h"
#include "Graphene.h"

/* MoS2 class */

class MoS2 : public Material
{
  public:
    // Constructor
    MoS2(const int rank, Parameters &rParameters, const Graphene &rGraphene):
      mRank(rank),
      mHoppingFile(rParameters.mHoppingFile)
      {
        mName = "MoS2";
        mNumWannier = 11;
        mNumWannierMo = 5;
        mNumNearestNeighbors = 7;
        mNumUnitCellX = 2;
        mNumUnitCellY = 3; 

        // Determines MoS2 lattice constant from graphene to exactly match the
        // widths of their supercells
        mLatticeConstant = rGraphene.GetLatticeConstant()
          * double(rGraphene.GetNumUnitCellY()) / double(mNumUnitCellY);
        mpHoppingData = (double *)calloc((2 * mNumNearestNeighbors + 1) * (2
              * mNumNearestNeighbors + 1) * mNumWannier * mNumWannier,
            sizeof(double));

        SetFermiLevel(rParameters);

        SetHoppingData();
      }

    ~MoS2(); // Destructor

//    // Sets the order of nearest neighbors
//    void SetNumNearestNeighbors(const int n_nearest_neighbors);
 
    // Gets the order of nearest neighbors
    const int GetNumNearestNeighbors() const;

    // Gets the hopping parameters data
    double *GetHoppingData() const;

    // Gets the Hamiltonian within MoS2 principal layer
    virtual cx_mat GetPrincipalLayerHamiltonian(const double phase) const;

    // Sets the interaction with the next MoS2 principal layer along the x axis
    virtual cx_mat GetPrincipalLayerInteractionX(const double phase) const;

  private:

    // Rank
    const int mRank;

    // Number of Wannier functions for Mo atom
    int mNumWannierMo;

    // Name of the Wannier90 output file containing MoS2 hopping parameters
    const std::string mHoppingFile;
    
    // MoS2 order of nearest neighbors in the x-y plane included in the
    // Wannier90 output file
    int mNumNearestNeighbors;

    // MoS2 hopping parameters
    double *mpHoppingData;

    // Sets MoS2 hopping parameters from the Wannier90 output file
    void SetHoppingData();

    // Sets MoS2 Fermi level according to input doping level
    void SetFermiLevel(Parameters &rParameters);

    // Used for constructing Hamiltonian between MoS2 unit cells
    cx_mat ConstructMoS2Hamiltonian(const int x, const int y, const double
        phase) const;
};

//inline 
//void MoS2::SetNumNearestNeighbors(const int n_nearest_neighbors)
//{
//  mNumNearestNeighbors = n_nearest_neighbors;
//}
  
inline
const int MoS2::GetNumNearestNeighbors() const
{
  return mNumNearestNeighbors;
}

inline
double *MoS2::GetHoppingData() const
{
  return mpHoppingData;
}

#endif
