// @author Wushi Dong

// PoissonSolver.cpp

#include "PoissonSolver.h"


PoissonSolver::~PoissonSolver(){}


void PoissonSolver::Run(PotentialProfile potential_profile, const ChargeProfile
    charge_profile, const FixedChargeProfile fixed_charge_profile, bool
    &rIsInnerConv)
{
  // Initialization
  double *pot = potential_profile.GetPotential();
  double *charge = charge_profile.GetCharge();
  double *charge_fixed = fixed_charge_profile.GetCharge();
  int n_grid = charge_profile.GetNumGrid();
  int n_grid_principal_layer = 6;

  // Potential grid spacing
  double DELTA_X_Gr = mLatticeConstantGraphene / sqrt(3.0) / 2.0;
  double DELTA_X_MoS2 = mLatticeConstantMoS2 / sqrt(3.0) / 2.0;
  
  // Jacobian
  mat J(n_grid, n_grid); 
  vec F(n_grid);
  // Potential changes
  vec det_pot(n_grid); 

  // Local judge for possion solver convergence
  int judge = 0;
  // Count of iteration runs
  int count = 0;
//  double sign;

  // Records initial potential
  double pot_initial[n_grid];
  for(int i = 0; i < n_grid; ++i)
    pot_initial[i] = pot[i];

//  if(mRank == 0)
//  {
//    cout << "Charge distribution (A^-1):" << endl;
//    for(int i = 0; i < n_grid; ++i)
//    {
//      cout << charge_fixed[i] - charge[i] << "  ";
//    }
//    cout << endl;
//  }

  while(judge == 0)
  {
    judge = 1;
    count++;

    if(count > mPoissonRunMax)
      break;

    J.zeros();
    F.zeros();
    det_pot.zeros();

    if(mRank == 0)
    {
      cout << endl << "Inner cycle number "<< count << ":" << endl;
    }

//    if(mRank == 0)
//    {
//      cout << "Updated charge prediction:" << endl;
//      for(int i = 0; i < n_grid; ++i)
//        cout << i << '\t' << charge[i] * exp((pot[i] - pot_initial[i])) <<
//        endl;;
//      cout << endl;
//    }

    // Updates F (Neumman boundary condition) Graphene
    F(0) = DELTA_X_Gr * DIL_Graphene * pot[1] - 2.0 * DELTA_X_Gr * DIL_Graphene
      * pot[0] + DELTA_X_Gr * DIL_Graphene * pot[1] - DELTA_X_Gr * DELTA_X_Gr
      * (charge[0] * exp((pot[0] - pot_initial[0]) / mKT) - charge_fixed[0]);
    for(int i = 1; i < mNumPrincipalLayerGraphene * n_grid_principal_layer - 2;
        ++i) { F(i) = DELTA_X_Gr * DIL_Graphene * pot[i - 1] - 2.0 * DELTA_X_Gr
      * DIL_Graphene * pot[i] + DELTA_X_Gr * DIL_Graphene * pot[i + 1]
        - DELTA_X_Gr * DELTA_X_Gr * (charge[i] * exp((pot[i] - pot_initial[i])
              / mKT) - charge_fixed[i]); }

    // The last graphene carbon slice
    F(mNumPrincipalLayerGraphene * n_grid_principal_layer - 2) = 0.5
      * mDistanceGrS2 * DIL_Graphene * pot[mNumPrincipalLayerGraphene
      * n_grid_principal_layer - 3] - (0.5 * mDistanceGrS2 * DIL_Graphene
          + DELTA_X_Gr * DIL_Graphene) * pot[mNumPrincipalLayerGraphene
      * n_grid_principal_layer - 2] + DELTA_X_Gr * DIL_Graphene
      * pot[mNumPrincipalLayerGraphene * n_grid_principal_layer - 1]
      - DELTA_X_Gr * 0.5 * mDistanceGrS2 * (charge[mNumPrincipalLayerGraphene
          * n_grid_principal_layer - 2] * exp((pot[mNumPrincipalLayerGraphene
              * n_grid_principal_layer - 2]
              - pot_initial[mNumPrincipalLayerGraphene * n_grid_principal_layer
              - 2]) / mKT) - charge_fixed[mNumPrincipalLayerGraphene
          * n_grid_principal_layer - 2]);
    
    // Gr - MoS2 interface
    F(mNumPrincipalLayerGraphene * n_grid_principal_layer - 1) = 0.5
      * mDistanceGrS2 * DIL_Graphene * pot[mNumPrincipalLayerGraphene
      * n_grid_principal_layer - 2] - (0.5 * mDistanceGrS2 * DIL_Graphene + 0.5
          * mDistanceGrS2 * DIL_MoS2) * pot[mNumPrincipalLayerGraphene
      * n_grid_principal_layer - 1] + 0.5 * mDistanceGrS2 * DIL_MoS2
      * pot[mNumPrincipalLayerGraphene * n_grid_principal_layer];

    // S2
    F(mNumPrincipalLayerGraphene * n_grid_principal_layer) = DELTA_X_MoS2
      * DIL_MoS2 * pot[mNumPrincipalLayerGraphene * n_grid_principal_layer - 1]
      - (DELTA_X_MoS2 * DIL_MoS2 + 0.5 * mDistanceGrS2 * DIL_MoS2)
      * pot[mNumPrincipalLayerGraphene * n_grid_principal_layer] + 0.5
      * mDistanceGrS2 * DIL_MoS2 * pot[mNumPrincipalLayerGraphene
      * n_grid_principal_layer + 1] - 0.5 * mDistanceGrS2 * DELTA_X_MoS2
      * (charge[mNumPrincipalLayerGraphene * n_grid_principal_layer]
          * exp((pot[mNumPrincipalLayerGraphene * n_grid_principal_layer]
              - pot_initial[mNumPrincipalLayerGraphene
              * n_grid_principal_layer]) / mKT)
          - charge_fixed[mNumPrincipalLayerGraphene
          * n_grid_principal_layer]);

    // MoS2
    for(int i = mNumPrincipalLayerGraphene * n_grid_principal_layer + 1;
        i < n_grid - 1; ++i) { F(i) = DELTA_X_MoS2 * DIL_MoS2 * pot[i - 1] - 2.0
      * DELTA_X_MoS2 * DIL_MoS2 * pot[i] + DELTA_X_MoS2 * DIL_MoS2 * pot[i + 1]
        - DELTA_X_MoS2 * DELTA_X_MoS2 * (charge[i] * exp((pot[i]
                - pot_initial[i]) / mKT) - charge_fixed[i]); } F(n_grid - 1)
        = DELTA_X_MoS2 * DIL_MoS2 * pot[n_grid - 2] - 2.0 * DELTA_X_MoS2
        * DIL_MoS2 * pot[n_grid - 1] + DELTA_X_MoS2 * DIL_MoS2 * pot[n_grid - 2]
        - DELTA_X_MoS2 * DELTA_X_MoS2 * (charge[n_grid - 1] * exp((pot[n_grid
                - 1] - pot_initial[n_grid - 1]) / mKT) - charge_fixed[n_grid
            - 1]);

//    if(mRank == 0)
//    {
//      cout << "F:" << endl;
//      for(int i = 0; i < n_grid; ++i)
//        cout << i << '\t' << F(i) << endl;
//      cout << endl;
//    }

    // Updates Jacobian
    // Graphene
    J(0, 0) = -2.0 * DELTA_X_Gr * DIL_Graphene;
    J(0, 1) = DELTA_X_Gr * DIL_Graphene + DELTA_X_Gr * DIL_Graphene; 
    for(int i = 1; i < mNumPrincipalLayerGraphene * n_grid_principal_layer - 2;
        ++i)
    {
      J(i, i) = -2.0 * DELTA_X_Gr * DIL_Graphene - DELTA_X_Gr * DELTA_X_Gr
        * charge[i - 1] * exp((pot[i] - pot_initial[i]) / mKT) / mKT; 
      J(i, i - 1) = DELTA_X_Gr * DIL_Graphene; 
      J(i, i + 1) = DELTA_X_Gr * DIL_Graphene; 
    }

    // The last graphene carbon slice
    J(mNumPrincipalLayerGraphene * n_grid_principal_layer - 2,
        mNumPrincipalLayerGraphene * n_grid_principal_layer - 2) = -(0.5
          * mDistanceGrS2 * DIL_Graphene + DELTA_X_Gr * DIL_Graphene) - DELTA_X_Gr
          * 0.5 * mDistanceGrS2 * charge[mNumPrincipalLayerGraphene
          * n_grid_principal_layer - 2] * exp((pot[mNumPrincipalLayerGraphene
                * n_grid_principal_layer - 2]
                - pot_initial[mNumPrincipalLayerGraphene
                * n_grid_principal_layer - 2]) / mKT) / mKT;
    J(mNumPrincipalLayerGraphene * n_grid_principal_layer - 2,
        mNumPrincipalLayerGraphene * n_grid_principal_layer - 3) = 0.5
      * mDistanceGrS2 * DIL_Graphene; 
    J(mNumPrincipalLayerGraphene * n_grid_principal_layer - 2,
        mNumPrincipalLayerGraphene * n_grid_principal_layer - 1) = DELTA_X_Gr
      * DIL_Graphene;

    // Gr - MoS2 interface
    J(mNumPrincipalLayerGraphene * n_grid_principal_layer - 1,
        mNumPrincipalLayerGraphene * n_grid_principal_layer - 1) = -(0.5
          * mDistanceGrS2 * DIL_Graphene + 0.5 * mDistanceGrS2 * DIL_MoS2);
    J(mNumPrincipalLayerGraphene * n_grid_principal_layer - 1,
        mNumPrincipalLayerGraphene * n_grid_principal_layer - 2) = 0.5
      * mDistanceGrS2 * DIL_Graphene; 
    J(mNumPrincipalLayerGraphene * n_grid_principal_layer - 1,
        mNumPrincipalLayerGraphene * n_grid_principal_layer) = 0.5
      * mDistanceGrS2 * DIL_MoS2;

    // S2
    J(mNumPrincipalLayerGraphene * n_grid_principal_layer,
        mNumPrincipalLayerGraphene * n_grid_principal_layer) = -(DELTA_X_MoS2
          * DIL_MoS2 + 0.5 * mDistanceGrS2 * DIL_MoS2) - 0.5 * mDistanceGrS2
          * DELTA_X_MoS2 * charge[mNumPrincipalLayerGraphene
          * n_grid_principal_layer] * exp((pot[mNumPrincipalLayerGraphene
                * n_grid_principal_layer]
                - pot_initial[mNumPrincipalLayerGraphene
                * n_grid_principal_layer]) / mKT) / mKT;
    J(mNumPrincipalLayerGraphene * n_grid_principal_layer,
        mNumPrincipalLayerGraphene * n_grid_principal_layer - 1) = DELTA_X_MoS2
      * DIL_MoS2; 
    J(mNumPrincipalLayerGraphene * n_grid_principal_layer,
        mNumPrincipalLayerGraphene * n_grid_principal_layer + 1) = 0.5
      * mDistanceGrS2 * DIL_MoS2;

    // MoS2
    for(int i = mNumPrincipalLayerGraphene * n_grid_principal_layer + 1;
        i < n_grid - 1; ++i) 
    {
      J(i, i) = -2.0 * DELTA_X_MoS2 * DIL_MoS2 - DELTA_X_MoS2 * DELTA_X_MoS2
        * charge[i] * exp((pot[i] - pot_initial[i]) / mKT) / mKT; 
      J(i, i - 1) = DELTA_X_MoS2 * DIL_MoS2; 
      J(i, i + 1) = DELTA_X_MoS2 * DIL_MoS2; 
    } 
    J(n_grid - 1, n_grid - 1) = -2.0 * DELTA_X_MoS2 * DIL_MoS2 - DELTA_X_MoS2
      * DELTA_X_MoS2 * charge[n_grid - 1] * exp((pot[n_grid - 1]
            - pot_initial[n_grid - 1]) / mKT) / mKT;
    J(n_grid - 1, n_grid - 2) = DELTA_X_MoS2 * DIL_MoS2 + DELTA_X_MoS2
      * DIL_MoS2; /* Neumman boundary condition */

//    if(mRank == 0)
//    {
//      cout << "J:" << endl;
//      cout << 0 << '\t' << J(0, 0) << '\t' << J(0, 1) << endl;
//      for(int i = 1; i < n_grid - 1; ++i)
//        cout << i << '\t' << J(i, i - 1) << '\t' << J(i, i) << '\t' << J(i, i + 1) << endl;
//      cout << n_grid - 1 << '\t' << J(n_grid - 1,n_grid - 2) << '\t' << J(n_grid - 1, n_grid - 1) << endl;
//      cout << endl;
//    }

    // Gets potential change
    det_pot = -J.i() * F;

//    if(mRank == 0)
//    {
//      cout << "Potential updates:" << endl;
//      for(int i = 0; i < n_grid; ++i)
//        cout << i << '\t' << det_pot(i) << endl;
//      cout << endl;
//    }

//    // Brown and Lindsay fix if necessary
//    for(int i = 0; i < n_grid; ++i)
//    {
//      if(abs(det_pot(i)) > 1)
//      {
//        if(mRank == 0)        
//          cout << "Fixed!" << endl;
//        if(det_pot(i) > 0)
//          sign = 1.0;
//        else
//          sign = -1.0;
//        if(abs(det_pot(i)) < 3.7)
//          det_pot(i) = sign * pow(abs(det_pot(i)), 0.2);
//        else
//          det_pot(i) = sign * log(abs(det_pot(i)));
//      }
//    }

    // Updates potential with damping
    for(int i = 0; i < n_grid; ++i)
    {
      pot[i] += (1.0 - mPotentialDampingPoisson) * det_pot(i);
    }
//    if(mRank == 0)
//    {    
//      cout << "Updated potential at current iteration (V):" << endl;
//      for (int i = 0; i < n_grid; ++i)
//        cout << i << "\t" << setprecision(7) << pot[i] << endl;
//    }

    // Judges convergence
    if(mRank == 0)
    {
      cout << "Max det_pot = " << (abs(det_pot)).max() << endl;
    }
    if(abs(det_pot).max() > mPotentialConv)
      judge = 0;
  }

  if(count > mPoissonRunMax)
  {
    if(mRank == 0)
    {
      cout << endl << "... Inner Poisson iteration NOT converged!" << endl;
    }
    rIsInnerConv = false;
  }
  else
  {
    if(mRank == 0)
    {    
      cout << "... converged at inner iteration number: " << count << endl << endl;
      cout << "Converged potential from Poisson solver (V):" << endl;
      for (int i = 0; i < n_grid; ++i)
        cout << i << "\t" << std::setprecision(7) << pot[i] << endl;
    }
  }

}


//// Testing
//#include "ChargeSolver.h"
//int main()
//{
//  // Input parameters
//  Parameters parameters;
//  parameters.ParseInputFile();
//  parameters.Print();  
//
//  DeviceEdgeContact device_edge_contact(parameters);
//
//  ChargeProfile charge_profile(30);
//  PotentialProfile potential_profile(30);
//  FixedChargeProfile fixed_charge(parameters);
//  fixed_charge.SetCharge(charge_profile.GetCharge());
//  double charge[30];
//  charge[30] = 0.0001;
//  charge_profile.SetCharge(charge);
//
//  cout << "At the beginning ..." << endl;
//  potential_profile.Print();
//  charge_profile.Print();
//  cout << endl;
//
//  cout << "Attaching the PoissonSolver ...  " << endl;
//  PoissonSolver poisson_solver(parameters, device_edge_contact);
//  cout << "Running the PoissonSolver ...  " << endl;
//  bool IsInnerConv = true;
//  bool rIsInnerConv = &IsInnerConv;
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

