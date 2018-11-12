
// @author Wushi Dong

// DeviceEdgeContact.h 
// Purpose: The electron transport device based on the edge contact. It consists
// of two leads, a central region and the materials that produce these parts.
// This is the main input object to the transport simulator.


#ifndef DEVICEEDGECONTACT_H
#define DEVICEEDGECONTACT_H

#include "Parameters.h"
#include "Material.h"
#include "Graphene.h"
#include "MoS2.h"
#include "Lead.h"
#include "CentralRegionEdgeContact.h"

class DeviceEdgeContact
{
  public:
    DeviceEdgeContact(const Parameters &rParameters, const Graphene &rGraphene,
        const MoS2 &rMoS2):
      mrGraphene(rGraphene),
      mLeftGrapheneLead(rParameters, rGraphene),
      mrMoS2(rMoS2),
      mRightMoS2Lead(rParameters, rMoS2),
      mCentralRegionEdgeContact(rParameters, rGraphene, rMoS2){
      }

    ~DeviceEdgeContact();

    // Graphene lead
    // Gets the Number of unit cell repetitions in a supercell along the x axis
    // for Graphene
    int GetNumUnitCellXGraphene() const;

    // Gets the Number of unit cell repetitions in a supercell along the Y axis
    // for Graphene
    int GetNumUnitCellYGraphene() const;
    
    // Gets the principal layer Hamiltonian within the left Graphene lead
    cx_mat GetPrincipalLayerHamiltonianGrapheneLead(const double phase) const;

    // Gets the principal layer interaction within the left Graphene lead
    cx_mat GetPrincipalLayerInteractionGrapheneLead(const double phase) const;
    
//    // Gets the interaction with the left Graphene lead
//    cx_mat GetInteractionWithLeftGrapheneLead(const double phase) const;

    // Gets the number of WFs in a single unit cell for Graphene
    const int GetNumWannierGraphene() const;

    // Gets the number of principal layers of the left Graphene region
    const int GetNumPrincipalLayerGraphene() const;

    // Gets the Graphene lattice constant
    double GetGrapheneLatticeConstant() const;

    
    // Central region
    // Gets the interaction with the left Graphene region
    cx_mat GetInteractionGrapheneS2(const double phase) const;

    // Gets the Hamiltonian of the interface S2 dimer
    cx_mat GetHamiltonianS2(const double phase) const;

    // Gets the number of WFs in a single unit cell for S2
    const int GetNumWannierS2() const;
    
    // Gets the number of WFs in a single unit cell for Mo
    const int GetNumWannierMo() const;
      
    // Gets the interaction with the right MoS2 region
    cx_mat GetInteractionS2MoS2(const double phase) const;


    // MoS2 lead
    // Gets the Number of unit cell repetitions in a supercell along the x axis
    // for MoS2
    int GetNumUnitCellXMoS2() const;

    // Gets the Number of unit cell repetitions in a supercell along the y axis
    // for MoS2
    int GetNumUnitCellYMoS2() const;

    // Gets the number of WFs in a single unit cell for MoS2
    const int GetNumWannierMoS2() const;

    // Gets the number of principal layers of the left Graphene region
    const int GetNumPrincipalLayerMoS2() const;

//    // Gets the interaction with the right MoS2 lead
//    cx_mat GetInteractionWithRightMoS2Lead(const double phase) const;

    // Gets the principal layer Hamiltonian within the right MoS2 lead
    cx_mat GetPrincipalLayerHamiltonianMoS2Lead(const double phase) const;

    // Gets the principal layer interaction within the right MoS2 lead
    cx_mat GetPrincipalLayerInteractionMoS2Lead(const double phase) const;

    // Gets the MoS2 lattice constant
    double GetMoS2LatticeConstant() const;

  private:

    const Graphene &mrGraphene;

    const Lead mLeftGrapheneLead;

    const MoS2 &mrMoS2;
    
    const Lead mRightMoS2Lead;

    const CentralRegionEdgeContact mCentralRegionEdgeContact;

};

inline
int DeviceEdgeContact::GetNumUnitCellXGraphene() const
{
  return mrGraphene.GetNumUnitCellX();
}

inline
int DeviceEdgeContact::GetNumUnitCellYGraphene() const
{
  return mrGraphene.GetNumUnitCellY();
}

inline 
cx_mat DeviceEdgeContact::GetPrincipalLayerHamiltonianGrapheneLead(const double
    phase) const
{
  return mLeftGrapheneLead.GetPrincipalLayerHamiltonian(phase);
}

inline 
cx_mat DeviceEdgeContact::GetPrincipalLayerInteractionGrapheneLead(const double
    phase) const
{
  return mLeftGrapheneLead.GetPrincipalLayerInteractionX(phase);
}

inline 
const int DeviceEdgeContact::GetNumWannierGraphene() const
{
  return mrGraphene.GetNumWannier();
}

inline
const int DeviceEdgeContact::GetNumPrincipalLayerGraphene() const
{
  return mCentralRegionEdgeContact.GetNumPrincipalLayerGraphene();
}

inline 
double DeviceEdgeContact::GetGrapheneLatticeConstant() const
{
  return mCentralRegionEdgeContact.GetLatticeConstantGraphene();
}


inline
cx_mat DeviceEdgeContact::GetInteractionGrapheneS2(double phase) const
{
  return mCentralRegionEdgeContact.GetInteractionGrapheneS2(phase);
}

inline
cx_mat DeviceEdgeContact::GetHamiltonianS2(double phase) const
{
  return mCentralRegionEdgeContact.GetHamiltonianS2(phase);
}

inline
const int DeviceEdgeContact::GetNumWannierS2() const
{
  return mCentralRegionEdgeContact.GetNumWannierS2();
}

inline
const int DeviceEdgeContact::GetNumWannierMo() const
{
  return mCentralRegionEdgeContact.GetNumWannierMo();
}

inline
cx_mat DeviceEdgeContact::GetInteractionS2MoS2(double phase) const
{
  return mCentralRegionEdgeContact.GetInteractionS2MoS2(phase);
}


inline
int DeviceEdgeContact::GetNumUnitCellXMoS2() const
{
  return mrMoS2.GetNumUnitCellX();
}

inline
int DeviceEdgeContact::GetNumUnitCellYMoS2() const
{
  return mrMoS2.GetNumUnitCellY();
}

inline 
const int DeviceEdgeContact::GetNumWannierMoS2() const
{
  return mrMoS2.GetNumWannier();
}

inline
const int DeviceEdgeContact::GetNumPrincipalLayerMoS2() const
{
  return mCentralRegionEdgeContact.GetNumPrincipalLayerMoS2();
}

inline cx_mat DeviceEdgeContact::GetPrincipalLayerHamiltonianMoS2Lead(const
    double phase) const
{
  return mRightMoS2Lead.GetPrincipalLayerHamiltonian(phase);
}

inline cx_mat DeviceEdgeContact::GetPrincipalLayerInteractionMoS2Lead(const
    double phase) const
{
  return mRightMoS2Lead.GetPrincipalLayerInteractionX(phase);
}

inline 
double DeviceEdgeContact::GetMoS2LatticeConstant() const
{
  return mCentralRegionEdgeContact.GetLatticeConstantMoS2();
}

//inline 
//cx_mat DeviceEdgeContact::GetInteractionWithLeftGrapheneLead(const double phase)
//  const
//{
//  return mLeftGrapheneLead.GetPrincipalLayerInteractionX(phase);
//}

//inline 
//cx_mat DeviceEdgeContact::GetInteractionWithRightMoS2Lead(const double phase)
//  const
//{
//  return mRightMoS2Lead.GetPrincipalLayerInteractionX(phase);
//}

#endif
