// @author Wushi Dong

// Graphene.h
// Purpose: Class that contains information for zigzag graphene

#ifndef GRAPHENE_H
#define GRAPHENE_H

#include "header.h"
#include "Parameters.h"
#include "Material.h"

// Contains information (name, geometry, Hamiltonian) of zigzag graphene
class Graphene : public Material
{
  public:
    // Constructor
    // @param rParameters Reference to a Parameters object
    Graphene(const Parameters &);

    // Destructor
    ~Graphene();

    // Sets the Graphene nearest hopping
    void SetHopping(const double hopping);

    // Sets the Hamiltonian within Graphene unit cell
    void SetUnitCellHamiltonian();

    // Sets the interaction with the next Graphene unit cell along the x axis
    void SetUnitCellInteractionX();

    // Sets the interaction with the next Graphene unit cell along the y axis
    void SetUnitCellInteractionY();

    // Gets the Hamiltonian within Graphene unit cell
    cx_mat GetUnitCellHamiltonian() const;

    // Gets the interaction with the next Graphene unit cell along the x axis
    cx_mat GetUnitCellInteractionX() const;

    // Gets the interaction with the next Graphene unit cell along the y axis
    cx_mat GetUnitCellInteractionY() const;

    // Gets the Hamiltonian within Graphene principal layer
    virtual cx_mat GetPrincipalLayerHamiltonian(const double phase) const;

    // Gets the interaction with the next Graphene principal layer along the
    // x axis
    virtual cx_mat GetPrincipalLayerInteractionX(const double phase) const;
    
  private:
    // Graphene nearest hopping
    double mHopping;

    // Hamiltonian within Material unit cell
    cx_mat mUnitCellHamiltonian;

    // Interaction with the next Material unit cell along the x axis
    cx_mat mUnitCellInteractionX;

    // Interaction with the next Material unit cell along the y axis
    cx_mat mUnitCellInteractionY;
};


inline
void Graphene::SetHopping(const double hopping)
{
  mHopping = hopping;
}

inline
cx_mat Graphene::GetUnitCellHamiltonian() const
{
  return mUnitCellHamiltonian;
}
inline
cx_mat Graphene::GetUnitCellInteractionX() const
{
  return mUnitCellInteractionX;
}

inline
cx_mat Graphene::GetUnitCellInteractionY() const
{
  return mUnitCellInteractionY;
}


#endif
