// @author Wushi Dong

// Material.h 
// Purpose: Class that contains information for a generic material

#ifndef MATERIAL_H 
#define MATERIAL_H

#include <string>
#include <armadillo>

using namespace arma;

// Contains information (name, geometry, Hamiltonian) of a generic material
class Material 
{ 
  public:
    // Constructor
    Material();

    // Destructor
    virtual ~Material();

    // Gets the material name
    const std::string GetName() const;

    // Gets the lattice constant
    const double GetLatticeConstant() const;

    // Gets the number of Wannier Funtions (WFs) used in a single unit cell
    const int GetNumWannier() const;

    // Gets the number of unit cell repetitions in a supercell along the x axis
    const int GetNumUnitCellX() const;

    // Gets the number of unit cell repetitions in a supercell along the x axis
    const int GetNumUnitCellY() const;

    // Gets the Hamiltonian within a principal layer
    //
    // @param phase The phase in the periodic direction for the current
    // calculation.  
    // @return Hamiltonian for a principal layer
    virtual cx_mat GetPrincipalLayerHamiltonian(const double phase) const
      = 0;

    // Gets the interaction between two adjacent principal layers
    //
    // @param phase The phase in the periodic direction for the current
    // calculation.
    // @return The interaction between two adjacent principal layers
    virtual cx_mat GetPrincipalLayerInteractionX(const double phase) const
      = 0;

  protected:
    // Material prefix
    std::string mName;

    // Lattice constant
    double mLatticeConstant;

    // Fermi level
    double mFermiLevel;

    // Number of WFs in a single unit cell
    int mNumWannier;

    // Number of unit cell repetitions in a supercell along the x axis
    int mNumUnitCellX;

    // Number of unit cell repetitions in a supercell along the y axis
    int mNumUnitCellY;
    
//    // Hamiltonian for one supercell 
//    cx_mat mPrincipalLayerHamiltonian;

//    // Interaction between two supercells along x-axis cx_mat
//    mPrincipalLayerInteractionX;

};


inline const std::string Material::GetName() const 
{ 
  return mName; 
}

inline const double Material::GetLatticeConstant() const 
{
  return mLatticeConstant; 
}

inline const int Material::GetNumWannier() const 
{
  return mNumWannier; 
}

inline const int Material::GetNumUnitCellX() const 
{
  return mNumUnitCellX; 
}

inline const int Material::GetNumUnitCellY() const 
{
  return mNumUnitCellY; 
}


#endif
