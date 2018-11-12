// @author Wushi Dong

// Lead.h 
// Purpose: Lead used by a electron transport device

#ifndef LEAD_H
#define LEAD_H

#include "Parameters.h"
#include "Material.h"

class Lead
{
  public:
    // Constructor
    Lead(const Parameters, const Material &rMaterial):
      mrMaterial(rMaterial){}; 
    ~Lead(); // Destructor

//    // Set Lead Hamiltonian
//    void SetLeadHamiltonian(const cx_mat LeadHamiltonian);

//    // Set Lead Interaction
//    void SetLeadInteraction(const cx_mat LeadInteraction);

//    // Set convergence criteria for calculating surface Green's functions from
//    input parameters
//    void SetCriteria(const Parameters parameters);

    // Get Lead Hamiltonian for one supercell
    arma::cx_mat GetPrincipalLayerHamiltonian(const double phase) const;

    // Get Lead Interaction between supercells
    arma::cx_mat GetPrincipalLayerInteractionX(const double phase) const;

//    // Calculate Lead surface Green's function by Sancho-Rubio's method
//    cx_mat GetSurfaceGreensFunction(const complex<double> energy, const double
//    eta, const double phase) const;
 
  private:

    // Material for making the lead
    const Material &mrMaterial;

//    // Lead Hamiltonian for one supercell
//    cx_mat mLeadHamiltonian;

//    // Lead Interaction between supercells
//    cx_mat mLeadInteraction;
   
};


//inline
//void Lead::SetLeadHamiltonian(const cx_mat LeadHamiltonian)
//{
//  mLeadHamiltonian = LeadHamiltonian;
//}
//
//inline
//void Lead::SetLeadInteraction(const cx_mat LeadInteraction)
//{
//  mLeadInteraction = LeadInteraction;
//}

//inline
//void Lead::SetCriteria(const Parameters parameters)
//{
//  mCriteria = parameters.mSurfaceGreensFunctionConv;
//}

inline
cx_mat Lead::GetPrincipalLayerHamiltonian(const double phase) const
{
  return mrMaterial.GetPrincipalLayerHamiltonian(phase);;
}

inline
cx_mat Lead::GetPrincipalLayerInteractionX(const double phase) const
{
  return mrMaterial.GetPrincipalLayerInteractionX(phase);;
}  

#endif
