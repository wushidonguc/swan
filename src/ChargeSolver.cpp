// @author Wushi Dong

// ChargeSolver.cpp

#include "ChargeSolver.h"


ChargeSolver::~ChargeSolver(){}


void ChargeSolver::SetHamiltonians(const double phase)
{
  mPrincipalLayerHamiltonianGraphene
    = mrDeviceEdgeContact.GetPrincipalLayerHamiltonianGrapheneLead(phase);

  mPrincipalLayerInteractionGraphene
    = mrDeviceEdgeContact.GetPrincipalLayerInteractionGrapheneLead(phase);

  //       = mrDeviceEdgeContact.GetInteractionWithLeftGrapheneLead(phase));

  mInteractionGrapheneS2 = mrDeviceEdgeContact.GetInteractionGrapheneS2(phase);

  mHamiltonianS2 = mrDeviceEdgeContact.GetHamiltonianS2(phase);

  mInteractionS2MoS2 = mrDeviceEdgeContact.GetInteractionS2MoS2(phase);

  //     = mrDeviceEdgeContact.GetInteractionWithRightMoS2Lead(phase));

  mPrincipalLayerHamiltonianMoS2
    = mrDeviceEdgeContact.GetPrincipalLayerHamiltonianMoS2Lead(phase);

  mPrincipalLayerInteractionMoS2
    = mrDeviceEdgeContact.GetPrincipalLayerInteractionMoS2Lead(phase);

}


cx_mat ChargeSolver::add_pot_Gr_left(cx_mat PL, double LEFT_POT)
{
  cx_mat PL_pot(PL);
  for(int i = 0; i < mNumWannierSupercellGraphene; ++i)
    PL_pot(i, i) -= LEFT_POT;
  return PL_pot;
}


cx_mat ChargeSolver::add_pot_Gr(cx_mat PL, int i_scat_Gr, double * pot)
{
  cx_mat PL_pot(PL);
  for(int i_rep = 0; i_rep < mNumUnitCellYGraphene; ++i_rep)
  {
    PL_pot(0 + i_rep * mNumWannierGraphene * mNumUnitCellXGraphene, 0 + i_rep
        * mNumWannierGraphene * mNumUnitCellXGraphene) -= pot[0 + i_scat_Gr
      * mNumGridPrincipalLayer]; 
    PL_pot(1 + i_rep * mNumWannierGraphene * mNumUnitCellXGraphene, 1 + i_rep
        * mNumWannierGraphene * mNumUnitCellXGraphene) -= pot[1 + i_scat_Gr
      * mNumGridPrincipalLayer]; 
    PL_pot(2 + i_rep * mNumWannierGraphene * mNumUnitCellXGraphene, 2 + i_rep
        * mNumWannierGraphene * mNumUnitCellXGraphene) -= pot[3 + i_scat_Gr
      * mNumGridPrincipalLayer]; 
    PL_pot(3 + i_rep * mNumWannierGraphene * mNumUnitCellXGraphene, 3 + i_rep
        * mNumWannierGraphene * mNumUnitCellXGraphene) -= pot[4 + i_scat_Gr
      * mNumGridPrincipalLayer];
  }
  return PL_pot;
}


cx_mat ChargeSolver::add_pot_s2(cx_mat PL, double * pot, int
    mNumPrincipalLayerGraphene)
{
  cx_mat PL_pot(PL);
  for(int i = 0; i < mNumWannierSupercellS2; ++i)
    PL_pot(i, i) -= pot[mNumPrincipalLayerGraphene * mNumGridPrincipalLayer];
  return PL_pot;
}


cx_mat ChargeSolver::add_pot_MoS2(cx_mat PL, int i_scat_MoS2, double * pot, int
    mNumPrincipalLayerGraphene)
{
  cx_mat PL_pot(PL);
  for(int i_rep = 0; i_rep < mNumUnitCellYMoS2; ++i_rep)
  {
    for(int i = 0; i < mNumWannierS2 - 1; ++i)
      PL_pot(i + i_rep * mNumWannierMoS2 * mNumUnitCellXMoS2, i + i_rep
          * mNumWannierMoS2 * mNumUnitCellXMoS2) -= pot[0 + i_scat_MoS2
        * mNumGridPrincipalLayer + mNumPrincipalLayerGraphene
        * mNumGridPrincipalLayer + 1];
    for(int i = mNumWannierS2 - 1; i < mNumWannierMoS2; ++i)
      PL_pot(i + i_rep * mNumWannierMoS2 * mNumUnitCellXMoS2, i + i_rep
          * mNumWannierMoS2 * mNumUnitCellXMoS2) -= pot[1 + i_scat_MoS2
        * mNumGridPrincipalLayer + mNumPrincipalLayerGraphene
        * mNumGridPrincipalLayer + 1];
    for(int i = mNumWannierMoS2; i < mNumWannierMoS2  + mNumWannierS2 - 1; ++i)
      PL_pot(i + i_rep * mNumWannierMoS2 * mNumUnitCellXMoS2, i + i_rep
          * mNumWannierMoS2 * mNumUnitCellXMoS2) -= pot[3 + i_scat_MoS2
        * mNumGridPrincipalLayer + mNumPrincipalLayerGraphene
        * mNumGridPrincipalLayer + 1];
    for(int i = mNumWannierMoS2  + mNumWannierS2 - 1; i < mNumWannierMoS2
        * mNumUnitCellXMoS2; ++i)
      PL_pot(i + i_rep * mNumWannierMoS2 * mNumUnitCellXMoS2, i + i_rep
          * mNumWannierMoS2 * mNumUnitCellXMoS2) -= pot[4 + i_scat_MoS2
        * mNumGridPrincipalLayer + mNumPrincipalLayerGraphene
        * mNumGridPrincipalLayer + 1];
  }
  return PL_pot;
}


cx_mat ChargeSolver::add_pot_MoS2_right(cx_mat PL, double RIGHT_POT)
{
  cx_mat PL_pot(PL);
  for(int i = 0; i < mNumWannierSupercellMoS2; ++i)
    PL_pot(i, i) -= RIGHT_POT;
  return PL_pot;
}


// Uses Sancho-Rubio's decimation method [Sancho84, Sancho85]
cx_mat ChargeSolver::GetLeadSurfaceGreensFunc(cx_mat
    principal_layer_hamiltonian, cx_mat principal_layer_interaction,
    std::complex<double> energy)
{
  cx_mat es = principal_layer_hamiltonian;
  cx_mat ee = principal_layer_hamiltonian;
  cx_mat a = principal_layer_interaction;
  cx_mat b = principal_layer_interaction.t();

  int nRows = es.n_rows;
  int nCols = es.n_cols;
  cx_mat E = (energy + mEta * constants::I) * eye(nRows, nCols);
 
  cx_mat es_temp(nRows, nCols);
  cx_mat ee_temp(nRows, nCols);
  cx_mat a_temp(nRows, nCols);
  cx_mat b_temp(nRows, nCols);  
  cx_mat g00(nRows, nCols);

  E = (energy + mEta * constants::I) * eye(nRows, nCols);

  while ((abs(a)).max() > mCriteria or (abs(b)).max() > mCriteria)
  {
    a_temp = a, b_temp = b, es_temp = es, ee_temp = ee;
    a = a_temp * (E - ee_temp).i() * a_temp;
    b = b_temp * (E - ee_temp).i() * b_temp;
    es = es_temp + a_temp * (E - ee_temp).i() * b_temp;
    ee = ee_temp + a_temp * (E - ee_temp).i() * b_temp + b_temp * (E
        - ee_temp).i() * a_temp;
  }

  g00 = (E - es).i();

  return g00;
}


void ChargeSolver::SolveCharge(ChargeProfile charge_profile, const
    PotentialProfile potential_profile)
{
  // Initialization
  double *charge = charge_profile.GetCharge();
  double *pot = potential_profile.GetPotential();
  int n_grid = potential_profile.GetNumGrid();

  // Surface Green's function of the left lead
  cx_mat g001(mNumWannierSupercellGraphene, mNumWannierSupercellGraphene);
  // Surface Green's function of the right lead
  cx_mat g002(mNumWannierSupercellMoS2, mNumWannierSupercellMoS2);
  // Retarded self-energy of the left lead
  cx_mat sigma1(mNumWannierSupercellGraphene, mNumWannierSupercellGraphene);
  // Retarded self-energy of the right lead
  cx_mat sigma2(mNumWannierSupercellMoS2, mNumWannierSupercellMoS2); 
  // Hamiltonian of the central region
  cx_mat H(mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
      + mNumWannierSupercellS2 + mNumPrincipalLayerMoS2
      * mNumWannierSupercellMoS2, mNumPrincipalLayerGraphene
      * mNumWannierSupercellGraphene + mNumWannierSupercellS2
      + mNumPrincipalLayerMoS2 * mNumWannierSupercellMoS2);
  // Retarded Green's function of the central region
  cx_mat G_r(size(H));
  // Retarded self-energy of the left lead (same size of the central region)
  cx_mat sigma1_full(size(H));
  // Retarded self-energy of the right lead (same size of the central region)
  cx_mat sigma2_full(size(H));
  // Broadening function of the left lead (same size of the central region)
  cx_mat gamma1_full(size(H)); 
  // Broadening function of the right lead (same size of the central region)
  cx_mat gamma2_full(size(H)); 

  // Diagonal terms of the retarded Green's function of graphene in the central
  // region
  cx_cube G_r_nn_Gr(mNumWannierSupercellGraphene, mNumWannierSupercellGraphene,
      mNumPrincipalLayerGraphene); 
  G_r_nn_Gr.zeros();

  // Diagonal terms of the retarded Green's function of S2 in the central region
  cx_mat G_r_nn_s2(mNumWannierSupercellS2, mNumWannierSupercellS2); 
  G_r_nn_s2.zeros();

  // Diagonal terms of the retarded Green's function of MoS2 in the central
  // region
  cx_cube G_r_nn_MoS2(mNumWannierSupercellMoS2, mNumWannierSupercellMoS2, mNumPrincipalLayerMoS2); 
  G_r_nn_MoS2.zeros();

  // Angle of current point on the contour integral
  double theta;
  // Number of sample points on the contour
  int N_theta = 5e1;
  // Angle step
  double step_theta = M_PI / (double)N_theta;
  // Complex energy of specific point on the contour
  std::complex<double> energy;
  // Lower bound of the contour (must be lower than the band bottom)
  double E_min = mEnergyRange.mStart - 1.0;
  // Lower electrochemical potential
  double mu_min = -10.0 * mKT;
  // Higher electrochemical potential
  double mu_max = +10.0 * mKT;
  // Center of the contour
  double center = 0.5 * (mu_min + E_min);
  // Radius of the contour
  double radius = 0.5 * (mu_min - E_min);
  // Value of the contour intefral
  double int_contour[n_grid];
  for(int i = 0; i < n_grid; ++i)
    int_contour[i] = 0.0;
  // Index of lower energy point
  int n_min = int((mu_min - mEnergyRange.mStart) / mEnergyRange.mSpacing);
  // Index of higher energy point
  int n_max = int((mu_max - mEnergyRange.mStart) / mEnergyRange.mSpacing);
  // Left and right spectral functions ...
  cx_mat A1(size(H));
  cx_mat A2(size(H));
  // And their real part
  double a1, a2;
  
  
  /* Equilibrium part */

  // Loops over semi-circle path
  for(int i_theta = 0; i_theta < N_theta; ++i_theta)
  {
    theta = step_theta * (double(i_theta) + 0.5);
    energy = center - radius * cos(theta) + constants::I * radius * sin(theta);

    // Surface Green's functions
    g001.zeros(); 
    g001
      = GetLeadSurfaceGreensFunc((add_pot_Gr_left(mPrincipalLayerHamiltonianGraphene,
              pot[0])).t(), mPrincipalLayerInteractionGraphene.t(), energy);
    g002.zeros(); 
    g002
      = GetLeadSurfaceGreensFunc(add_pot_MoS2_right(mPrincipalLayerHamiltonianMoS2,
            pot[n_grid - 1]), mPrincipalLayerInteractionMoS2, energy);

    // Retarded self-energies from contacts
    sigma1.zeros(); sigma1 = mPrincipalLayerInteractionGraphene.t() * g001
      * mPrincipalLayerInteractionGraphene; sigma2.zeros(); 
    sigma2 = mPrincipalLayerInteractionMoS2 * g002
      * mPrincipalLayerInteractionMoS2.t();

    // Next uses the recursive Green's function method 
    // Forward sweep 
    // Graphene layer
    G_r_nn_Gr.slice(0) = (energy * eye(mNumWannierSupercellGraphene,
          mNumWannierSupercellGraphene)
        - add_pot_Gr(mPrincipalLayerHamiltonianGraphene, 0, pot) - sigma1).i();
 
    for(int i_scat = 1; i_scat < mNumPrincipalLayerGraphene; ++i_scat)
    {
      G_r_nn_Gr.slice(i_scat) = (energy * eye(mNumWannierSupercellGraphene,
            mNumWannierSupercellGraphene)
          - add_pot_Gr(mPrincipalLayerHamiltonianGraphene, i_scat, pot)
          - mPrincipalLayerInteractionGraphene.t() * G_r_nn_Gr.slice(i_scat - 1)
          * mPrincipalLayerInteractionGraphene).i();
    }

    // S2 layer
    G_r_nn_s2 = (energy * eye(mNumWannierSupercellS2, mNumWannierSupercellS2)
        - add_pot_s2(mHamiltonianS2, pot, mNumPrincipalLayerGraphene)
        - mInteractionGrapheneS2.t()
        * G_r_nn_Gr.slice(mNumPrincipalLayerGraphene - 1)
        * mInteractionGrapheneS2).i();

    // MoS2 layer
    G_r_nn_MoS2.slice(0) = (energy * eye(mNumWannierSupercellMoS2,
          mNumWannierSupercellMoS2)
        - add_pot_MoS2(mPrincipalLayerHamiltonianMoS2, 0, pot,
          mNumPrincipalLayerGraphene) - mInteractionS2MoS2.t() * G_r_nn_s2
        * mInteractionS2MoS2).i();
    for(int i_scat = 1; i_scat < mNumPrincipalLayerMoS2 - 1; ++i_scat)
    {
      G_r_nn_MoS2.slice(i_scat) = (energy * eye(mNumWannierSupercellMoS2,
            mNumWannierSupercellMoS2)
          - add_pot_MoS2(mPrincipalLayerHamiltonianMoS2, i_scat, pot,
            mNumPrincipalLayerGraphene) - mPrincipalLayerInteractionMoS2.t()
          * G_r_nn_MoS2.slice(i_scat - 1) * mPrincipalLayerInteractionMoS2).i();
    }
    G_r_nn_MoS2.slice(mNumPrincipalLayerMoS2 - 1) = (energy
        * eye(mNumWannierSupercellMoS2, mNumWannierSupercellMoS2)
        - add_pot_MoS2(mPrincipalLayerHamiltonianMoS2, mNumPrincipalLayerMoS2
          - 1, pot, mNumPrincipalLayerGraphene)
        - mPrincipalLayerInteractionMoS2.t()
        * G_r_nn_MoS2.slice(mNumPrincipalLayerMoS2 - 2)
        * mPrincipalLayerInteractionMoS2 - sigma2).i();

    // Backward sweep
    // MoS2 layer
    for(int i_scat = mNumPrincipalLayerMoS2 - 2; i_scat > -1; --i_scat)
      G_r_nn_MoS2.slice(i_scat) = G_r_nn_MoS2.slice(i_scat)
        + G_r_nn_MoS2.slice(i_scat) * mPrincipalLayerInteractionMoS2
        * G_r_nn_MoS2.slice(i_scat + 1) * mPrincipalLayerInteractionMoS2.t()
        * G_r_nn_MoS2.slice(i_scat);

    // S2 layer
    G_r_nn_s2 = G_r_nn_s2 + G_r_nn_s2 * mInteractionS2MoS2
      * G_r_nn_MoS2.slice(0) * mInteractionS2MoS2.t() * G_r_nn_s2;

    // Graphene layer
    G_r_nn_Gr.slice(mNumPrincipalLayerGraphene - 1)
      = G_r_nn_Gr.slice(mNumPrincipalLayerGraphene - 1)
      + G_r_nn_Gr.slice(mNumPrincipalLayerGraphene - 1) * mInteractionGrapheneS2
      * G_r_nn_s2 * mInteractionGrapheneS2.t()
      * G_r_nn_Gr.slice(mNumPrincipalLayerGraphene - 1);
    for(int i_scat = mNumPrincipalLayerGraphene - 2; i_scat > -1; --i_scat)
      G_r_nn_Gr.slice(i_scat) = G_r_nn_Gr.slice(i_scat)
        + G_r_nn_Gr.slice(i_scat) * mPrincipalLayerInteractionGraphene
        * G_r_nn_Gr.slice(i_scat + 1) * mPrincipalLayerInteractionGraphene.t()
        * G_r_nn_Gr.slice(i_scat);


    // Contour integral

    // Graphene
    for(int i_scat = 0; i_scat < mNumPrincipalLayerGraphene; ++i_scat)
    {
      for(int j = 0; j < mNumWannierSupercellGraphene; ++j)
      {
        switch(j % 4)
        {
          case 0: 
            int_contour[0 + i_scat * mNumGridPrincipalLayer] +=
              imag(trace(G_r_nn_Gr.slice(i_scat).submat(j, j, j, j))
                  * radius * (sin(theta) + constants::I * cos(theta)) * step_theta); 
            break;
          case 1: 
            int_contour[1 + i_scat * mNumGridPrincipalLayer] +=
              imag(trace(G_r_nn_Gr.slice(i_scat).submat(j, j, j, j))
                  * radius * (sin(theta) + constants::I * cos(theta)) * step_theta); 
            break;
          case 2: 
            int_contour[3 + i_scat * mNumGridPrincipalLayer] +=
              imag(trace(G_r_nn_Gr.slice(i_scat).submat(j, j, j, j))
                  * radius * (sin(theta) + constants::I * cos(theta)) * step_theta); 
            break;
          case 3: 
            int_contour[4 + i_scat * mNumGridPrincipalLayer] +=
              imag(trace(G_r_nn_Gr.slice(i_scat).submat(j, j, j, j))
                  * radius * (sin(theta) + constants::I * cos(theta)) * step_theta); 
            break;
        }
      }
    }

    // Gr - MoS2 interface
    int_contour[mNumPrincipalLayerGraphene * mNumGridPrincipalLayer - 1] = 0.0;

    // S2
    int_contour[mNumPrincipalLayerGraphene * mNumGridPrincipalLayer] +=
      imag(trace(G_r_nn_s2) * radius * (sin(theta) + constants::I * cos(theta)) 
          * step_theta);

    // MoS2
    for(int i_scat = 0; i_scat < mNumPrincipalLayerMoS2; ++i_scat)
    {
      for(int i_rep = 0; i_rep < 3; ++i_rep)
      {
        int_contour[0 + i_scat * mNumGridPrincipalLayer
          + mNumPrincipalLayerGraphene * mNumGridPrincipalLayer + 1] +=
          imag(trace(G_r_nn_MoS2.slice(i_scat).submat(0 + i_rep
                  * mNumWannierMoS2 * mNumUnitCellXMoS2, 0 + i_rep
                  * mNumWannierMoS2 * mNumUnitCellXMoS2, mNumWannierS2
                  - 1 + i_rep * mNumWannierMoS2 * mNumUnitCellXMoS2 - 1,
                  mNumWannierS2 - 1 + i_rep * mNumWannierMoS2
                  * mNumUnitCellXMoS2 - 1)) * radius * (sin(theta) + 
                  constants::I * cos(theta)) * step_theta);
        int_contour[1 + i_scat * mNumGridPrincipalLayer
          + mNumPrincipalLayerGraphene * mNumGridPrincipalLayer + 1] +=
          imag(trace(G_r_nn_MoS2.slice(i_scat).submat(mNumWannierS2 - 1 + i_rep
                  * mNumWannierMoS2 * mNumUnitCellXMoS2, mNumWannierS2
                  - 1 + i_rep * mNumWannierMoS2 * mNumUnitCellXMoS2,
                  mNumWannierMoS2 + i_rep * mNumWannierMoS2 * mNumUnitCellXMoS2
                  - 1, mNumWannierMoS2 + i_rep * mNumWannierMoS2
                  * mNumUnitCellXMoS2 - 1)) * radius * (sin(theta) + 
                  constants::I * cos(theta)) * step_theta);
        int_contour[3 + i_scat * mNumGridPrincipalLayer
          + mNumPrincipalLayerGraphene * mNumGridPrincipalLayer + 1] +=
          imag(trace(G_r_nn_MoS2.slice(i_scat).submat(mNumWannierMoS2 + i_rep
                  * mNumWannierMoS2 * mNumUnitCellXMoS2, mNumWannierMoS2 + i_rep
                  * mNumWannierMoS2 * mNumUnitCellXMoS2, mNumWannierMoS2
                  + mNumWannierS2 - 1 + i_rep * mNumWannierMoS2
                  * mNumUnitCellXMoS2 - 1, mNumWannierMoS2 + mNumWannierS2
                  - 1 + i_rep * mNumWannierMoS2 * mNumUnitCellXMoS2 - 1))
              * radius * (sin(theta) + constants::I * cos(theta)) * step_theta);
        int_contour[4 + i_scat * mNumGridPrincipalLayer
          + mNumPrincipalLayerGraphene * mNumGridPrincipalLayer + 1] +=
          imag(trace(G_r_nn_MoS2.slice(i_scat).submat(mNumWannierMoS2
                  + mNumWannierS2 - 1 + i_rep * mNumWannierMoS2
                  * mNumUnitCellXMoS2, mNumWannierMoS2 + mNumWannierS2
                  - 1 + i_rep * mNumWannierMoS2 * mNumUnitCellXMoS2,
                  mNumWannierMoS2 * mNumUnitCellXMoS2 + i_rep * mNumWannierMoS2
                  * mNumUnitCellXMoS2 - 1, mNumWannierMoS2 * mNumUnitCellXMoS2
                  + i_rep * mNumWannierMoS2 * mNumUnitCellXMoS2 - 1)) * radius
              * (sin(theta) + constants::I * cos(theta)) * step_theta);
      }
    }
  }
  for(int i_Grid = 0; i_Grid < n_grid; ++i_Grid)
  {
    charge[i_Grid] += -2.0 / M_PI * int_contour[i_Grid];
  }

  
  /* Non-equilibrium part */

  // Resets Hamiltonian to be safe
  H.zeros();

  // Graphene layer
  H.submat(0, 0, mNumWannierSupercellGraphene - 1, mNumWannierSupercellGraphene
      - 1) = add_pot_Gr(mPrincipalLayerHamiltonianGraphene, 0, pot);
  for(int i_scat = 1; i_scat < mNumPrincipalLayerGraphene; ++i_scat)
  {
    H.submat(i_scat * mNumWannierSupercellGraphene, i_scat
        * mNumWannierSupercellGraphene, (i_scat + 1)
        * mNumWannierSupercellGraphene - 1, (i_scat + 1)
        * mNumWannierSupercellGraphene - 1)
      = add_pot_Gr(mPrincipalLayerHamiltonianGraphene, i_scat, pot);
    H.submat((i_scat - 1) * mNumWannierSupercellGraphene, i_scat
        * mNumWannierSupercellGraphene, i_scat * mNumWannierSupercellGraphene
        - 1, (i_scat + 1) * mNumWannierSupercellGraphene - 1)
      = mPrincipalLayerInteractionGraphene;
    H.submat(i_scat * mNumWannierSupercellGraphene, (i_scat - 1)
        * mNumWannierSupercellGraphene, (i_scat + 1)
        * mNumWannierSupercellGraphene - 1, i_scat
        * mNumWannierSupercellGraphene - 1)
      = mPrincipalLayerInteractionGraphene.t();
  }

  // S dimer layer
  H.submat(mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene,
      mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene,
      mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
      + mNumWannierSupercellS2 - 1, mNumPrincipalLayerGraphene
      * mNumWannierSupercellGraphene + mNumWannierSupercellS2 - 1)
    = add_pot_s2(mHamiltonianS2, pot, mNumPrincipalLayerGraphene);
  H.submat((mNumPrincipalLayerGraphene - 1) * mNumWannierSupercellGraphene,
      mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene,
      mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene - 1,
      mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
      + mNumWannierSupercellS2 - 1) = mInteractionGrapheneS2;
  H.submat(mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene,
      (mNumPrincipalLayerGraphene - 1) * mNumWannierSupercellGraphene,
      mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
      + mNumWannierSupercellS2 - 1, mNumPrincipalLayerGraphene
      * mNumWannierSupercellGraphene - 1) = mInteractionGrapheneS2.t();

  // MoS2 layer
  H.submat(mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
      + mNumWannierSupercellS2, mNumPrincipalLayerGraphene
      * mNumWannierSupercellGraphene + mNumWannierSupercellS2,
      mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
      + mNumWannierSupercellS2 + mNumWannierSupercellMoS2 - 1,
      mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
      + mNumWannierSupercellS2 + mNumWannierSupercellMoS2 - 1)
    = add_pot_MoS2(mPrincipalLayerHamiltonianMoS2, 0, pot,
        mNumPrincipalLayerGraphene); 
  H.submat(mNumPrincipalLayerGraphene
          * mNumWannierSupercellGraphene, mNumPrincipalLayerGraphene
          * mNumWannierSupercellGraphene + mNumWannierSupercellS2,
          mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
          + mNumWannierSupercellS2 - 1, mNumPrincipalLayerGraphene
          * mNumWannierSupercellGraphene + mNumWannierSupercellS2
          + mNumWannierSupercellMoS2 - 1) = mInteractionS2MoS2;
  H.submat(mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
      + mNumWannierSupercellS2, mNumPrincipalLayerGraphene
      * mNumWannierSupercellGraphene, mNumPrincipalLayerGraphene
      * mNumWannierSupercellGraphene + mNumWannierSupercellS2
      + mNumWannierSupercellMoS2 - 1, mNumPrincipalLayerGraphene
      * mNumWannierSupercellGraphene + mNumWannierSupercellS2 - 1)
    = mInteractionS2MoS2.t();
  for(int i_scat = 1; i_scat < mNumPrincipalLayerMoS2; ++i_scat)
  {
    H.submat(mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
        + mNumWannierSupercellS2 + i_scat * mNumWannierSupercellMoS2,
        mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
        + mNumWannierSupercellS2 + i_scat * mNumWannierSupercellMoS2,
        mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
        + mNumWannierSupercellS2 + (i_scat + 1) * mNumWannierSupercellMoS2 - 1,
        mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
        + mNumWannierSupercellS2 + (i_scat + 1) * mNumWannierSupercellMoS2 - 1)
      = add_pot_MoS2(mPrincipalLayerHamiltonianMoS2, i_scat, pot,
          mNumPrincipalLayerGraphene); 
    H.submat(mNumPrincipalLayerGraphene
            * mNumWannierSupercellGraphene + mNumWannierSupercellS2 + (i_scat
              - 1) * mNumWannierSupercellMoS2, mNumPrincipalLayerGraphene
            * mNumWannierSupercellGraphene + mNumWannierSupercellS2 + i_scat
            * mNumWannierSupercellMoS2, mNumPrincipalLayerGraphene
            * mNumWannierSupercellGraphene + mNumWannierSupercellS2 + i_scat
            * mNumWannierSupercellMoS2 - 1, mNumPrincipalLayerGraphene
            * mNumWannierSupercellGraphene + mNumWannierSupercellS2 + (i_scat
              + 1) * mNumWannierSupercellMoS2 - 1)
            = mPrincipalLayerInteractionMoS2;
    H.submat(mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
        + mNumWannierSupercellS2 + i_scat * mNumWannierSupercellMoS2,
        mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
        + mNumWannierSupercellS2 + (i_scat - 1) * mNumWannierSupercellMoS2,
        mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
        + mNumWannierSupercellS2 + (i_scat + 1) * mNumWannierSupercellMoS2 - 1,
        mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
        + mNumWannierSupercellS2 + i_scat * mNumWannierSupercellMoS2 - 1)
      = mPrincipalLayerInteractionMoS2.t();
  }

  // Loops over energy grids
  for(int i = n_min; i < n_max; ++i)
  {
    double energy  = mEnergyRange.mStart + (double(i) + 0.5)
      * mEnergyRange.mSpacing;
    double f1 = mUtils.GetFermiDistribution(energy, mVoltageBias * 0.5, mKT);
    double f2 = mUtils.GetFermiDistribution(energy, -mVoltageBias * 0.5, mKT);

    // Surface Green's functions
    g001.zeros(); 
    g001
      = GetLeadSurfaceGreensFunc((add_pot_Gr_left(mPrincipalLayerHamiltonianGraphene,
              pot[0])).t(), mPrincipalLayerInteractionGraphene.t(), energy);
    g002.zeros(); 
    g002
      = GetLeadSurfaceGreensFunc(add_pot_MoS2_right(mPrincipalLayerHamiltonianMoS2,
            pot[n_grid - 1]), mPrincipalLayerInteractionMoS2, energy);

    // Retarded and advanced self-energies from contacts
    sigma1_full.zeros(); 
    sigma1_full.submat(0, 0, mNumWannierSupercellGraphene - 1,
        mNumWannierSupercellGraphene - 1)
      = mPrincipalLayerInteractionGraphene.t() * g001
      * mPrincipalLayerInteractionGraphene; 
    sigma2_full.zeros();
    sigma2_full.submat(mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
        + mNumWannierSupercellS2 + (mNumPrincipalLayerMoS2 - 1)
        * mNumWannierSupercellMoS2, mNumPrincipalLayerGraphene
        * mNumWannierSupercellGraphene + mNumWannierSupercellS2
        + (mNumPrincipalLayerMoS2 - 1) * mNumWannierSupercellMoS2,
        mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
        + mNumWannierSupercellS2 + mNumPrincipalLayerMoS2
        * mNumWannierSupercellMoS2 - 1, mNumPrincipalLayerGraphene
        * mNumWannierSupercellGraphene + mNumWannierSupercellS2
        + mNumPrincipalLayerMoS2 * mNumWannierSupercellMoS2 - 1)
      = mPrincipalLayerInteractionMoS2 * g002
      * mPrincipalLayerInteractionMoS2.t();

    // Direct inversion to get retarded Green's function
    G_r.zeros();
    G_r = (energy * eye(size(G_r)) - H - sigma1_full - sigma2_full).i();

    // Gamma from the contacts
    gamma1_full.zeros();
    gamma1_full = (sigma1_full - sigma1_full.t()) * constants::I;
    gamma2_full.zeros();
    gamma2_full = (sigma2_full - sigma2_full.t()) * constants::I;

    // Left and right spectral functions
    A1 = G_r * gamma1_full * G_r.t();
    A2 = G_r * gamma2_full * G_r.t();

    // Calculates electron density from left and right spectral funcitons
    // Gr
    for(int i_scat = 0; i_scat < mNumPrincipalLayerGraphene; ++i_scat)
    {
      for(int j = 0; j < mNumWannierSupercellGraphene; ++j)
      {
        switch(j % 4)
        {
          case 0:
            a1 = real(A1(i_scat * mNumGridPrincipalLayer + j, i_scat
                  * mNumGridPrincipalLayer + j)); 
            a2 = real(A2(i_scat * mNumGridPrincipalLayer + j, i_scat
                  * mNumGridPrincipalLayer + j)); 
            charge[0 + i_scat * mNumGridPrincipalLayer] += 1.0 / M_PI * (a1 * f1
                + a2 * f2) * mEnergyRange.mSpacing;
            break;
          case 1:
            a1 = real(A1(i_scat * mNumGridPrincipalLayer + j, i_scat
                  * mNumGridPrincipalLayer + j)); 
            a2 = real(A2(i_scat * mNumGridPrincipalLayer + j, i_scat
                  * mNumGridPrincipalLayer + j)); 
            charge[1 + i_scat * mNumGridPrincipalLayer] += 1.0 / M_PI * (a1 * f1
                + a2 * f2) * mEnergyRange.mSpacing;
            break;
          case 2:
            a1 = real(A1(i_scat * mNumGridPrincipalLayer + j, i_scat
                  * mNumGridPrincipalLayer + j)); a2 = real(A2(i_scat
                    * mNumGridPrincipalLayer + j, i_scat
                    * mNumGridPrincipalLayer + j)); charge[3 + i_scat
                  * mNumGridPrincipalLayer] += 1.0 / M_PI * (a1 * f1 + a2 * f2)
                  * mEnergyRange.mSpacing;
            break;
          case 3:
            a1 = real(A1(i_scat * mNumGridPrincipalLayer + j, i_scat
                  * mNumGridPrincipalLayer + j)); a2 = real(A2(i_scat
                    * mNumGridPrincipalLayer + j, i_scat
                    * mNumGridPrincipalLayer + j)); charge[4 + i_scat
                  * mNumGridPrincipalLayer] += 1.0 / M_PI * (a1 * f1 + a2 * f2)
                  * mEnergyRange.mSpacing;
            break;
        }
      }
    }

    // Gr - MoS2 interface
    charge[mNumPrincipalLayerGraphene * mNumGridPrincipalLayer - 1] = 0.0;

    // S2
    a1 = real(trace(A1.submat(mNumPrincipalLayerGraphene
            * mNumWannierSupercellGraphene, mNumPrincipalLayerGraphene
            * mNumWannierSupercellGraphene, mNumPrincipalLayerGraphene
            * mNumWannierSupercellGraphene + mNumWannierSupercellS2 - 1,
            mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
            + mNumWannierSupercellS2 - 1))); 
    a2 = real(trace(A2.submat(mNumPrincipalLayerGraphene
              * mNumWannierSupercellGraphene, mNumPrincipalLayerGraphene
              * mNumWannierSupercellGraphene, mNumPrincipalLayerGraphene
              * mNumWannierSupercellGraphene + mNumWannierSupercellS2 - 1,
              mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
              + mNumWannierSupercellS2 - 1))); 
    charge[mNumPrincipalLayerGraphene * mNumGridPrincipalLayer] += 1.0 / M_PI
      * (a1 * f1 + a2 * f2) * mEnergyRange.mSpacing;

    // MoS2
    for(int i_scat = 0; i_scat < mNumPrincipalLayerMoS2; ++i_scat)
    {
      for(int i_rep = 0; i_rep < 3; ++i_rep)
      {
        a1 = real(trace(A1.submat(0 + i_rep * mNumWannierMoS2
                * mNumUnitCellXMoS2 + mNumPrincipalLayerGraphene
                * mNumWannierSupercellGraphene + mNumWannierSupercellS2 + i_scat
                * mNumWannierSupercellMoS2, 0 + i_rep * mNumWannierMoS2
                * mNumUnitCellXMoS2 + mNumPrincipalLayerGraphene
                * mNumWannierSupercellGraphene + mNumWannierSupercellS2 + i_scat
                * mNumWannierSupercellMoS2, mNumWannierS2 - 1 + i_rep
                * mNumWannierMoS2 * mNumUnitCellXMoS2
                + mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
                + mNumWannierSupercellS2 + i_scat * mNumWannierSupercellMoS2
                - 1, mNumWannierS2 - 1 + i_rep * mNumWannierMoS2
                * mNumUnitCellXMoS2 + mNumPrincipalLayerGraphene
                * mNumWannierSupercellGraphene + mNumWannierSupercellS2 + i_scat
                * mNumWannierSupercellMoS2 - 1)));
        a2 = real(trace(A2.submat(0
                  + i_rep * mNumWannierMoS2 * mNumUnitCellXMoS2
                  + mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
                  + mNumWannierSupercellS2 + i_scat * mNumWannierSupercellMoS2,
                  0 + i_rep * mNumWannierMoS2 * mNumUnitCellXMoS2
                  + mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
                  + mNumWannierSupercellS2 + i_scat * mNumWannierSupercellMoS2,
                  mNumWannierS2 - 1 + i_rep * mNumWannierMoS2
                  * mNumUnitCellXMoS2 + mNumPrincipalLayerGraphene
                  * mNumWannierSupercellGraphene + mNumWannierSupercellS2
                  + i_scat * mNumWannierSupercellMoS2 - 1, mNumWannierS2
                  - 1 + i_rep * mNumWannierMoS2 * mNumUnitCellXMoS2
                  + mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
                  + mNumWannierSupercellS2 + i_scat * mNumWannierSupercellMoS2
                  - 1))); 
        charge[0 + i_scat * mNumGridPrincipalLayer + mNumPrincipalLayerGraphene 
          * mNumGridPrincipalLayer + 1] += 1.0 / M_PI * (a1 * f1 + a2 * f2)
          * mEnergyRange.mSpacing;

        a1 = real(trace(A1.submat(mNumWannierS2 - 1 + i_rep * mNumWannierMoS2
                * mNumUnitCellXMoS2 + mNumPrincipalLayerGraphene
                * mNumWannierSupercellGraphene + mNumWannierSupercellS2 + i_scat
                * mNumWannierSupercellMoS2, mNumWannierS2 - 1 + i_rep
                * mNumWannierMoS2 * mNumUnitCellXMoS2
                + mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
                + mNumWannierSupercellS2 + i_scat * mNumWannierSupercellMoS2,
                mNumWannierMoS2 + i_rep * mNumWannierMoS2 * mNumUnitCellXMoS2
                + mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
                + mNumWannierSupercellS2 + i_scat * mNumWannierSupercellMoS2
                - 1, mNumWannierMoS2 + i_rep * mNumWannierMoS2
                * mNumUnitCellXMoS2 + mNumPrincipalLayerGraphene
                * mNumWannierSupercellGraphene + mNumWannierSupercellS2 + i_scat
                * mNumWannierSupercellMoS2 - 1))); 
        a2 = real(trace(A2.submat(mNumWannierS2 - 1 + i_rep * mNumWannierMoS2
                  * mNumUnitCellXMoS2 + mNumPrincipalLayerGraphene
                  * mNumWannierSupercellGraphene + mNumWannierSupercellS2
                  + i_scat * mNumWannierSupercellMoS2, mNumWannierS2 - 1 + i_rep
                  * mNumWannierMoS2 * mNumUnitCellXMoS2
                  + mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
                  + mNumWannierSupercellS2 + i_scat * mNumWannierSupercellMoS2,
                  mNumWannierMoS2 + i_rep * mNumWannierMoS2 * mNumUnitCellXMoS2
                  + mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
                  + mNumWannierSupercellS2 + i_scat * mNumWannierSupercellMoS2
                  - 1, mNumWannierMoS2 + i_rep * mNumWannierMoS2
                  * mNumUnitCellXMoS2 + mNumPrincipalLayerGraphene
                  * mNumWannierSupercellGraphene + mNumWannierSupercellS2
                  + i_scat * mNumWannierSupercellMoS2 - 1))); 
        charge[1 + i_scat * mNumGridPrincipalLayer + mNumPrincipalLayerGraphene
          * mNumGridPrincipalLayer + 1] += 1.0 / M_PI * (a1 * f1 + a2 * f2)
          * mEnergyRange.mSpacing;

        a1 = real(trace(A1.submat(mNumWannierMoS2 + i_rep * mNumWannierMoS2
                * mNumUnitCellXMoS2 + mNumPrincipalLayerGraphene
                * mNumWannierSupercellGraphene + mNumWannierSupercellS2 + i_scat
                * mNumWannierSupercellMoS2, mNumWannierMoS2 + i_rep
                * mNumWannierMoS2 * mNumUnitCellXMoS2
                + mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
                + mNumWannierSupercellS2 + i_scat * mNumWannierSupercellMoS2,
                mNumWannierMoS2  + mNumWannierS2 - 1 + i_rep * mNumWannierMoS2
                * mNumUnitCellXMoS2 + mNumPrincipalLayerGraphene
                * mNumWannierSupercellGraphene + mNumWannierSupercellS2 + i_scat
                * mNumWannierSupercellMoS2 - 1, mNumWannierMoS2  + mNumWannierS2
                - 1 + i_rep * mNumWannierMoS2 * mNumUnitCellXMoS2
                + mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
                + mNumWannierSupercellS2 + i_scat * mNumWannierSupercellMoS2
                - 1))); 
        a2 = real(trace(A2.submat(mNumWannierMoS2 + i_rep
                  * mNumWannierMoS2 * mNumUnitCellXMoS2
                  + mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
                  + mNumWannierSupercellS2 + i_scat * mNumWannierSupercellMoS2,
                  mNumWannierMoS2 + i_rep * mNumWannierMoS2 * mNumUnitCellXMoS2
                  + mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
                  + mNumWannierSupercellS2 + i_scat * mNumWannierSupercellMoS2,
                  mNumWannierMoS2  + mNumWannierS2 - 1 + i_rep * mNumWannierMoS2
                  * mNumUnitCellXMoS2 + mNumPrincipalLayerGraphene
                  * mNumWannierSupercellGraphene + mNumWannierSupercellS2
                  + i_scat * mNumWannierSupercellMoS2 - 1, mNumWannierMoS2
                  + mNumWannierS2 - 1 + i_rep * mNumWannierMoS2
                  * mNumUnitCellXMoS2 + mNumPrincipalLayerGraphene
                  * mNumWannierSupercellGraphene + mNumWannierSupercellS2
                  + i_scat * mNumWannierSupercellMoS2 - 1))); 
        charge[3 + i_scat * mNumGridPrincipalLayer + mNumPrincipalLayerGraphene
                * mNumGridPrincipalLayer + 1] += 1.0 / M_PI * (a1 * f1 + a2
                    * f2) * mEnergyRange.mSpacing;

        a1 = real(trace(A1.submat(mNumWannierMoS2  + mNumWannierS2 - 1 + i_rep
                * mNumWannierMoS2 * mNumUnitCellXMoS2
                + mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
                + mNumWannierSupercellS2 + i_scat * mNumWannierSupercellMoS2,
                mNumWannierMoS2  + mNumWannierS2 - 1 + i_rep * mNumWannierMoS2
                * mNumUnitCellXMoS2 + mNumPrincipalLayerGraphene
                * mNumWannierSupercellGraphene + mNumWannierSupercellS2 + i_scat
                * mNumWannierSupercellMoS2, mNumWannierMoS2 * mNumUnitCellXMoS2
                + i_rep * mNumWannierMoS2 * mNumUnitCellXMoS2
                + mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
                + mNumWannierSupercellS2 + i_scat * mNumWannierSupercellMoS2
                - 1, mNumWannierMoS2 * mNumUnitCellXMoS2 + i_rep
                * mNumWannierMoS2 * mNumUnitCellXMoS2
                + mNumPrincipalLayerGraphene * mNumWannierSupercellGraphene
                + mNumWannierSupercellS2 + i_scat * mNumWannierSupercellMoS2
                - 1))); 
        a2 = real(trace(A2.submat(mNumWannierMoS2
                  + mNumWannierS2 - 1 + i_rep * mNumWannierMoS2
                  * mNumUnitCellXMoS2 + mNumPrincipalLayerGraphene
                  * mNumWannierSupercellGraphene + mNumWannierSupercellS2
                  + i_scat * mNumWannierSupercellMoS2, mNumWannierMoS2
                  + mNumWannierS2 - 1 + i_rep * mNumWannierMoS2
                  * mNumUnitCellXMoS2 + mNumPrincipalLayerGraphene
                  * mNumWannierSupercellGraphene + mNumWannierSupercellS2
                  + i_scat * mNumWannierSupercellMoS2, mNumWannierMoS2
                  * mNumUnitCellXMoS2 + i_rep * mNumWannierMoS2
                  * mNumUnitCellXMoS2 + mNumPrincipalLayerGraphene
                  * mNumWannierSupercellGraphene + mNumWannierSupercellS2
                  + i_scat * mNumWannierSupercellMoS2 - 1, mNumWannierMoS2
                  * mNumUnitCellXMoS2 + i_rep * mNumWannierMoS2
                  * mNumUnitCellXMoS2 + mNumPrincipalLayerGraphene
                  * mNumWannierSupercellGraphene + mNumWannierSupercellS2
                  + i_scat * mNumWannierSupercellMoS2 - 1))); 
        charge[4 + i_scat * mNumGridPrincipalLayer + mNumPrincipalLayerGraphene
                * mNumGridPrincipalLayer + 1] += 1.0 / M_PI * (a1 * f1 + a2
                    * f2) * mEnergyRange.mSpacing;
      }
    }
  }
}


void ChargeSolver::SolveTransmissionAndLDOS(double * transmission, double * LDOS, const PotentialProfile &rConvergedPotentialProfile, const EnergyRange &rEnergyRangeDense)
{
  // Initialize
  double *pot = rConvergedPotentialProfile.GetPotential();
  int n_grid = rConvergedPotentialProfile.GetNumGrid();
  
  // Surface Green's function of the left lead
  cx_mat g001(mNumWannierSupercellGraphene, mNumWannierSupercellGraphene);
  // Surface Green's function of the right lead
  cx_mat g002(mNumWannierSupercellMoS2, mNumWannierSupercellMoS2);
  // Retarded self-energy of the left lead
  cx_mat sigma1(mNumWannierSupercellGraphene, mNumWannierSupercellGraphene);
  // Retarded self-energy of the right lead
  cx_mat sigma2(mNumWannierSupercellMoS2, mNumWannierSupercellMoS2); 
  // Broadening function of the left lead
  cx_mat gamma1(mNumWannierSupercellGraphene, mNumWannierSupercellGraphene);
  // Broadening function of the right lead
  cx_mat gamma2(mNumWannierSupercellMoS2, mNumWannierSupercellMoS2); 

  // Diagonal terms of the retarded Green's function of graphene in the central
  // region
  cx_cube G_r_nn_Gr(mNumWannierSupercellGraphene, mNumWannierSupercellGraphene,
      mNumPrincipalLayerGraphene); 
  G_r_nn_Gr.zeros();

  // Diagonal terms of the retarded Green's function of S2 in the central region
  cx_mat G_r_nn_s2(mNumWannierSupercellS2, mNumWannierSupercellS2); 
  G_r_nn_s2.zeros();

  // Diagonal terms of the retarded Green's function of MoS2 in the central
  // region
  cx_cube G_r_nn_MoS2(mNumWannierSupercellMoS2, mNumWannierSupercellMoS2, mNumPrincipalLayerMoS2);
  G_r_nn_MoS2.zeros();

  // Complex energy of specific point on the contour
  std::complex<double> energy;

  // Loops over energy grids
	for(int i = 0; i < rEnergyRangeDense.mN; ++i)
  {
		energy  = rEnergyRangeDense.mStart + (double(i) + 0.5) * rEnergyRangeDense.mSpacing;

    // Top-right block of the retarded Green's function matrix
    cx_mat G_0n;

    // Surface Green's functions
    g001.zeros(); 
    g001
      = GetLeadSurfaceGreensFunc((add_pot_Gr_left(mPrincipalLayerHamiltonianGraphene,
              pot[0])).t(), mPrincipalLayerInteractionGraphene.t(), energy);
    g002.zeros(); 
    g002
      = GetLeadSurfaceGreensFunc(add_pot_MoS2_right(mPrincipalLayerHamiltonianMoS2,
            pot[n_grid - 1]), mPrincipalLayerInteractionMoS2, energy);

    // Retarded self-energies from contacts
    sigma1.zeros();
    sigma1 = mPrincipalLayerInteractionGraphene.t() * g001
      * mPrincipalLayerInteractionGraphene;
    sigma2.zeros();
    sigma2 = mPrincipalLayerInteractionMoS2 * g002
      * mPrincipalLayerInteractionMoS2.t();

    // Broadening functions from contacts
    gamma1.zeros();
    gamma1 = (sigma1 - sigma1.t()) * constants::I;
    gamma2.zeros();
    gamma2 = (sigma2 - sigma2.t()) * constants::I;

    // Next uses the recursive Green's function method 
    // Forward sweep 
    // Graphene layer
    G_r_nn_Gr.slice(0) = (energy * eye(mNumWannierSupercellGraphene,
          mNumWannierSupercellGraphene)
        - add_pot_Gr(mPrincipalLayerHamiltonianGraphene, 0, pot) - sigma1).i();
    G_0n = G_r_nn_Gr.slice(0);
 
    for(int i_scat = 1; i_scat < mNumPrincipalLayerGraphene; ++i_scat)
    {
      G_r_nn_Gr.slice(i_scat) = (energy * eye(mNumWannierSupercellGraphene,
            mNumWannierSupercellGraphene)
          - add_pot_Gr(mPrincipalLayerHamiltonianGraphene, i_scat, pot)
          - mPrincipalLayerInteractionGraphene.t() * G_r_nn_Gr.slice(i_scat - 1)
          * mPrincipalLayerInteractionGraphene).i();
      G_0n = G_0n * mPrincipalLayerInteractionGraphene * G_r_nn_Gr.slice(i_scat);
    }

    // S2 layer
    G_r_nn_s2 = (energy * eye(mNumWannierSupercellS2, mNumWannierSupercellS2)
        - add_pot_s2(mHamiltonianS2, pot, mNumPrincipalLayerGraphene)
        - mInteractionGrapheneS2.t()
        * G_r_nn_Gr.slice(mNumPrincipalLayerGraphene - 1)
        * mInteractionGrapheneS2).i();
    G_0n = G_0n * mInteractionGrapheneS2 * G_r_nn_s2;

    // MoS2 layer
    G_r_nn_MoS2.slice(0) = (energy * eye(mNumWannierSupercellMoS2,
          mNumWannierSupercellMoS2)
        - add_pot_MoS2(mPrincipalLayerHamiltonianMoS2, 0, pot,
          mNumPrincipalLayerGraphene) - mInteractionS2MoS2.t() * G_r_nn_s2
        * mInteractionS2MoS2).i();
    G_0n = G_0n * mInteractionS2MoS2 * G_r_nn_MoS2.slice(0);
    for(int i_scat = 1; i_scat < mNumPrincipalLayerMoS2 - 1; ++i_scat)
    {
      G_r_nn_MoS2.slice(i_scat) = (energy * eye(mNumWannierSupercellMoS2,
            mNumWannierSupercellMoS2)
          - add_pot_MoS2(mPrincipalLayerHamiltonianMoS2, i_scat, pot,
            mNumPrincipalLayerGraphene) - mPrincipalLayerInteractionMoS2.t()
          * G_r_nn_MoS2.slice(i_scat - 1) * mPrincipalLayerInteractionMoS2).i();
      G_0n = G_0n * mPrincipalLayerInteractionMoS2 * G_r_nn_MoS2.slice(i_scat);
    }
    G_r_nn_MoS2.slice(mNumPrincipalLayerMoS2 - 1) = (energy
        * eye(mNumWannierSupercellMoS2, mNumWannierSupercellMoS2)
        - add_pot_MoS2(mPrincipalLayerHamiltonianMoS2, mNumPrincipalLayerMoS2
          - 1, pot, mNumPrincipalLayerGraphene)
        - mPrincipalLayerInteractionMoS2.t()
        * G_r_nn_MoS2.slice(mNumPrincipalLayerMoS2 - 2)
        * mPrincipalLayerInteractionMoS2 - sigma2).i();
    G_0n = G_0n * mPrincipalLayerInteractionMoS2 * G_r_nn_MoS2.slice(mNumPrincipalLayerMoS2 - 1);

    // Calculates transmission
    transmission[i] += (trace(gamma1 * G_0n * gamma2 * G_0n.t())).real();    

    // Backward sweep
    // MoS2 layer
    for(int i_scat = mNumPrincipalLayerMoS2 - 2; i_scat > -1; --i_scat)
      G_r_nn_MoS2.slice(i_scat) = G_r_nn_MoS2.slice(i_scat)
        + G_r_nn_MoS2.slice(i_scat) * mPrincipalLayerInteractionMoS2
        * G_r_nn_MoS2.slice(i_scat + 1) * mPrincipalLayerInteractionMoS2.t()
        * G_r_nn_MoS2.slice(i_scat);

    // S2 layer
    G_r_nn_s2 = G_r_nn_s2 + G_r_nn_s2 * mInteractionS2MoS2
      * G_r_nn_MoS2.slice(0) * mInteractionS2MoS2.t() * G_r_nn_s2;

    // Graphene layer
    G_r_nn_Gr.slice(mNumPrincipalLayerGraphene - 1)
      = G_r_nn_Gr.slice(mNumPrincipalLayerGraphene - 1)
      + G_r_nn_Gr.slice(mNumPrincipalLayerGraphene - 1) * mInteractionGrapheneS2
      * G_r_nn_s2 * mInteractionGrapheneS2.t()
      * G_r_nn_Gr.slice(mNumPrincipalLayerGraphene - 1);
    for(int i_scat = mNumPrincipalLayerGraphene - 2; i_scat > -1; --i_scat)
      G_r_nn_Gr.slice(i_scat) = G_r_nn_Gr.slice(i_scat)
        + G_r_nn_Gr.slice(i_scat) * mPrincipalLayerInteractionGraphene
        * G_r_nn_Gr.slice(i_scat + 1) * mPrincipalLayerInteractionGraphene.t()
        * G_r_nn_Gr.slice(i_scat);

    // calculate LDOS
    // Graphene
    for(int i_scat = 0; i_scat < mNumPrincipalLayerGraphene; ++i_scat)
    {
      for(int j = 0; j < 16; ++j)
      {
        switch(j % 4)
        {
          case 0: LDOS[(0 + i_scat * 4) * rEnergyRangeDense.mN + i] -= 2.0 / M_PI
                  * imag(trace(G_r_nn_Gr.slice(i_scat).submat(j, j, j, j)));
                  break;
          case 1: LDOS[(1 + i_scat * 4) * rEnergyRangeDense.mN + i] -= 2.0 / M_PI
                  * imag(trace(G_r_nn_Gr.slice(i_scat).submat(j, j, j, j)));
                  break;
          case 2: LDOS[(2 + i_scat * 4) * rEnergyRangeDense.mN + i] -= 2.0 / M_PI
                  * imag(trace(G_r_nn_Gr.slice(i_scat).submat(j, j, j, j)));
                  break;
          case 3: LDOS[(3 + i_scat * 4) * rEnergyRangeDense.mN + i] -= 2.0 / M_PI
                  * imag(trace(G_r_nn_Gr.slice(i_scat).submat(j, j, j, j)));
                  break;
        }
      }
    }

    // S dimer
    LDOS[(mNumPrincipalLayerGraphene * 4) * rEnergyRangeDense.mN + i] -= 2.0 / M_PI * imag(trace(G_r_nn_s2));

    // MoS2
    for(int i_scat = 0; i_scat < mNumPrincipalLayerMoS2; ++i_scat)
    {
      for(int i_rep = 0; i_rep < 3; ++i_rep)
      {
        LDOS[(0 + i_scat * 4 + mNumPrincipalLayerGraphene * 4 + 1) * rEnergyRangeDense.mN + i] -= 2.0 / M_PI
          * imag(trace(G_r_nn_MoS2.slice(i_scat).submat(0 + i_rep * 22,
                  0 + i_rep * 22, 5 + i_rep * 22 - 1, 5 + i_rep * 22 - 1)));
        LDOS[(1 + i_scat * 4 + mNumPrincipalLayerGraphene * 4 + 1) * rEnergyRangeDense.mN + i] -= 2.0 / M_PI
          * imag(trace(G_r_nn_MoS2.slice(i_scat).submat(5 + i_rep * 22,
                  5 + i_rep * 22, 11 + i_rep * 22 - 1, 11 + i_rep * 22 - 1)));
        LDOS[(2 + i_scat * 4 + mNumPrincipalLayerGraphene * 4 + 1) * rEnergyRangeDense.mN + i] -= 2.0 / M_PI
          * imag(trace(G_r_nn_MoS2.slice(i_scat).submat(11 + i_rep * 22, 11
                  + i_rep * 22, 16 + i_rep * 22 - 1, 16 + i_rep * 22 - 1)));
        LDOS[(3 + i_scat * 4 + mNumPrincipalLayerGraphene * 4 + 1) * rEnergyRangeDense.mN + i] -= 2.0 / M_PI
          * imag(trace(G_r_nn_MoS2.slice(i_scat).submat(16 + i_rep * 22, 16
                  + i_rep * 22, 22 + i_rep * 22 - 1, 22 + i_rep * 22 - 1)));
      }
    }  

  }
}


//// Testing
//#include "PoissonSolver.h"
//int main()
//{
//  // Initialization
////  PhaseRange phase_range(0., 0., 1);
//  double phase = 0.0;
//  
//  // Input parameters
//  Parameters parameters;
//  parameters.ParseInputFile();
//  parameters.Print();
//
//  DeviceEdgeContact device_edge_contact(parameters);
//
//  ChargeProfile charge_profile(parameters);
//  double *charge = charge_profile.GetCharge();
//  PotentialProfile potential_profile(parameters);
//  FixedChargeProfile fixed_charge(parameters);
//
//  cout << "At the beginning ..." << endl;
//  potential_profile.Print();
//  charge_profile.Print();
//  cout << endl;
//
//  cout << "Attaching the ChargeSolver ...  " << endl;
//  ChargeSolver charge_solver(parameters, device_edge_contact);
//  charge_solver.SetHamiltonians(phase);
//  cout << "Running the ChargeSolver ...  " << endl;
//  charge_solver.Run(charge, potential_profile);
//
//  cout << "After running ChargeSolver ..." << endl;
//  potential_profile.Print();
//  charge_profile.Print();
//  cout << endl;
//
//  cout << "Attaching the PoissonSolver ...  " << endl;
//  PoissonSolver poisson_solver(0, parameters, device_edge_contact);
//  cout << "Running the PoissonSolver ...  " << endl;
//  bool IsInnerConv = true;
//  bool *pIsInnerConv = &IsInnerConv;
//  poisson_solver.Run(potential_profile, charge_profile, fixed_charge, pIsInnerConv);
//
//  cout << "After running PoissonSolver ..." << endl;
//  cout << "PoissonSolver converged? " << IsInnerConv << endl;
//  potential_profile.Print();
//  charge_profile.Print();
//  cout << endl;
//  
//  return 0;
//}
