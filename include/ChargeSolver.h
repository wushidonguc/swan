// @author Wushi Dong

// ChargeSolver.h
// Purpose: Calculates charge profile based on the Keldysh formalism given
// certain electrostatic potential.

#ifndef CHARGERSOLVER_H
#define CHARGERSOLVER_H

#include "Parameters.h"
#include "DeviceEdgeContact.h"
#include "ChargeProfile.h"
#include "PotentialProfile.h"
#include "EnergyRange.h"
#include "Utils.h"

class ChargeSolver
{
  public:
    // Constructor    
    ChargeSolver(const Parameters &rParameters, const DeviceEdgeContact
        &rDeviceEdgeContact):

      mEnergyRange(rParameters),

      mrDeviceEdgeContact(rDeviceEdgeContact),

      mNumUnitCellXGraphene(rDeviceEdgeContact.GetNumUnitCellXGraphene()),

      mNumUnitCellYGraphene(rDeviceEdgeContact.GetNumUnitCellYGraphene()),

      mNumWannierGraphene(rDeviceEdgeContact.GetNumWannierGraphene()),

      mNumPrincipalLayerGraphene(
          rDeviceEdgeContact.GetNumPrincipalLayerGraphene()),

      mNumWannierS2(rDeviceEdgeContact.GetNumWannierS2()),

      mNumWannierMo(rDeviceEdgeContact.GetNumWannierMo()),

      mNumUnitCellXMoS2(rDeviceEdgeContact.GetNumUnitCellXMoS2()),

      mNumUnitCellYMoS2(rDeviceEdgeContact.GetNumUnitCellYMoS2()),

      mNumWannierMoS2(rDeviceEdgeContact.GetNumWannierMoS2()),

      mNumPrincipalLayerMoS2(rDeviceEdgeContact.GetNumPrincipalLayerMoS2()),
    
      mCriteria(rParameters.mSurfaceGreensFunctionConv),
      
      mKT(constants::BOLTZMANN_CONSTANT * rParameters.mTemperature),
      
      mVoltageBias(rParameters.mVoltageBias)

      {
        mNumGridPrincipalLayer = 6;
        mNumWannierSupercellGraphene = mNumWannierGraphene
          * mNumUnitCellXGraphene * mNumUnitCellYGraphene;
        mNumWannierSupercellS2 = mNumWannierS2 * 1 * mNumUnitCellYMoS2;
        mNumWannierSupercellMoS2 = mNumWannierMoS2 * mNumUnitCellXMoS2
          * mNumUnitCellYMoS2;
        mEta = mKT;
        mNumGridLDOS = 4;
      }

    // Destructor
    ~ChargeSolver();

    // Sets Hamiltonian from the member device object
    void SetHamiltonians(const double phase);

    // Solves for charge density profile
    void SolveCharge(ChargeProfile, const PotentialProfile);

    // Solves for transmission and LDOS
    void SolveTransmissionAndLDOS(double *, double *, const PotentialProfile &, const EnergyRange &);
    
  private:
    EnergyRange mEnergyRange;

    const DeviceEdgeContact &mrDeviceEdgeContact;

    Utils mUtils;
    
    int mNumGridPrincipalLayer;
    
    // Graphene
    int mNumUnitCellXGraphene;
    
    int mNumUnitCellYGraphene;

    int mNumWannierGraphene;

    int mNumWannierSupercellGraphene;

    int mNumPrincipalLayerGraphene;

    cx_mat mPrincipalLayerHamiltonianGraphene;

    cx_mat mPrincipalLayerInteractionGraphene;

    // Graphene - S2
    cx_mat mInteractionGrapheneS2;

    // S2
    int mNumWannierS2;

    int mNumWannierSupercellS2;

    int mNumWannierMo;

    cx_mat mHamiltonianS2;

    // S2 - MoS2
    cx_mat mInteractionS2MoS2;

    //MoS2
    int mNumUnitCellXMoS2;
    
    int mNumUnitCellYMoS2;

    int mNumWannierMoS2;
    
    int mNumWannierSupercellMoS2;

    int mNumPrincipalLayerMoS2;

    cx_mat mPrincipalLayerHamiltonianMoS2;

    cx_mat mPrincipalLayerInteractionMoS2;
    
    // Convergence criteria for calculating surface Green's functions
    double mCriteria;

    // Boltzman constant * Temperature
    double mKT;

    // Broadening
    double mEta;

    double mVoltageBias;

    // Transmission and LDOS calculation
    // Number of grids for LDOS
    int mNumGridLDOS;
    
    
    // Add electrostatic potential to the Hamiltonian of left graphene lead
    cx_mat add_pot_Gr_left(cx_mat PL, double LEFT_POT);
    
    // Add electrostatic potential to Hamiltonian of graphene in central region
    cx_mat add_pot_Gr(cx_mat PL, int i_scat_Gr, double * pot);

    // Add electrostatic potential to the Hamiltonian of interface S2
    cx_mat add_pot_s2(cx_mat PL, double * pot, int mNumPrincipalLayerGraphene);

    // Add electrostatic potential to Hamiltonian of MoS2 in central region
    cx_mat add_pot_MoS2(cx_mat PL, int i_scat_MoS2, double * pot, int
        mNumPrincipalLayerGraphene);

    // Add electrostatic potential to the Hamiltonian of right MoS2 lead
    cx_mat add_pot_MoS2_right(cx_mat PL, double RIGHT_POT);
    
    cx_mat GetLeadSurfaceGreensFunc(cx_mat principal_layer_hamiltonian, cx_mat principal_layer_interaction, std::complex<double> energy);

};


#endif
