
// @author Wushi Dong

// CentralRegionEdgeContact.h 
// Purpose: The central region of a edge contact device. It contains information
// on the geometry and Hamiltonian of the edge contact structure.


#ifndef CENTRALREGIONEDGECONTACT_H
#define CENTRALREGIONEDGECONTACT_H

#include "header.h"
#include "Parameters.h"
#include "Material.h"
#include "Graphene.h"
#include "MoS2.h"


class CentralRegionEdgeContact
{
  public:
    // Constructor
    CentralRegionEdgeContact(const Parameters &parameters, const Graphene
        &rGraphene, const MoS2 &rMoS2):
      mrGraphene(rGraphene),
      mNumPrincipalLayerGraphene(parameters.mNumPrincipalLayerGraphene),
      mNumUnitCellXGraphene(rGraphene.GetNumUnitCellX()),
      mNumUnitCellYGraphene(rGraphene.GetNumUnitCellY()),
      mNumWannierGraphene(rGraphene.GetNumWannier()),
      mLatticeConstantGraphene(rGraphene.GetLatticeConstant()),
      mrMoS2(rMoS2),
      mNumPrincipalLayerMoS2(parameters.mNumPrincipalLayerMoS2),
      mNumUnitCellXMoS2(rMoS2.GetNumUnitCellX()),
      mNumUnitCellYMoS2(rMoS2.GetNumUnitCellY()),
      mNumWannierMoS2(rMoS2.GetNumWannier()),
      mLatticeConstantMoS2(rMoS2.GetLatticeConstant()),
      mNumNearestNeighborsMoS2(rMoS2.GetNumNearestNeighbors()),
      mpHoppingData(rMoS2.GetHoppingData()),
      mYShift(0.2 * rGraphene.GetLatticeConstant()),  /* We find that such shift
                                                         can give the most
                                                         stable interfacial
                                                         structure */
      mHoppingCS(parameters.mHoppingCS),
      mDistanceGrS2(parameters.mDistanceGrS2)
      {
//        cout << "[Device] Creating central region for Graphene-MoS2 edge
//        contact ..." << endl << endl;

        mNumWannierS2 = 6;
        mNumWannierMo = 5;
      };

    // Destructor
    ~CentralRegionEdgeContact();

    // Gets the interaction with the left Graphene region
    cx_mat GetInteractionGrapheneS2(const double phase) const;

    // Gets the number of principal layers of the left Graphene region
    const int GetNumPrincipalLayerGraphene() const;

    // Gets the graphene lattice constant
    double GetLatticeConstantGraphene() const;

    // Gets the number of WFs in a unit cell of S2
    const int GetNumWannierS2() const;

    // Gets the number of WFs in a unit cell of Mo
    const int GetNumWannierMo() const;

    // Gets the Hamiltonian of the interface S2 dimer
    cx_mat GetHamiltonianS2(const double phase) const;

    // Gets the interaction with the right MoS2 region
    cx_mat GetInteractionS2MoS2(const double phase) const;

    // Gets the number of principal layers of the left Graphene region
    const int GetNumPrincipalLayerMoS2() const;

    // Gets the MoS2 lattice constant
    double GetLatticeConstantMoS2() const;

  private:
    const Graphene &mrGraphene;
    int mNumPrincipalLayerGraphene;
    int mNumUnitCellXGraphene;
    int mNumUnitCellYGraphene;
    int mNumWannierGraphene;
    int mLatticeConstantGraphene;

    // Total number of Wannier functions for the two S atom in MoS2 unit cell
    int mNumWannierS2;
    int mNumWannierMo;

    const MoS2 &mrMoS2;
    int mNumPrincipalLayerMoS2;
    int mNumUnitCellXMoS2;
    int mNumUnitCellYMoS2;
    int mNumWannierMoS2;
    int mLatticeConstantMoS2;
    int mNumNearestNeighborsMoS2;
    double *mpHoppingData;
    
    // Relative shift along the y-axis between the two materials
    double mYShift;

    // Hopping between the pz-like WF of the graphene carbon atom and the
    // px-like WF of the MoS2 sulfur atom closest to the interface
    double mHoppingCS;

    // Interfacial distance between graphene and MoS2
    double mDistanceGrS2;

    // Gets the interaction between interface S2 and the MoS2 unit cells
    cx_mat GetS2MoS2Interaction(const int y, const double phase) const;

};

inline
const int CentralRegionEdgeContact::GetNumPrincipalLayerGraphene() const
{
  return mNumPrincipalLayerGraphene;
}

inline
double CentralRegionEdgeContact::GetLatticeConstantGraphene() const
{
  return mrGraphene.GetLatticeConstant();
}

inline
const int CentralRegionEdgeContact::GetNumWannierS2() const
{
  return mNumWannierS2;
}

inline
const int CentralRegionEdgeContact::GetNumWannierMo() const
{
  return mNumWannierMo;
}

inline
const int CentralRegionEdgeContact::GetNumPrincipalLayerMoS2() const
{
  return mNumPrincipalLayerMoS2;
}

inline
double CentralRegionEdgeContact::GetLatticeConstantMoS2() const
{
  return mrMoS2.GetLatticeConstant();
}


#endif
