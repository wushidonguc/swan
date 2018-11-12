// @author Wushi Dong

// SelfConsistentSolver.cpp

#include "SelfConsistentSolver.h"


SelfConsistentSolver::~SelfConsistentSolver(){}


void SelfConsistentSolver::Run()
{
  // Initialize
  // Number of grid point for the charge and potential profile
  const int n_grid = mPotentialProfile.GetNumGrid();

  double *charge = mChargeProfile.GetCharge();

  double *charge_fixed = mFixedChargeProfile.GetCharge();

  double *charge_sum = (double *)calloc(n_grid, sizeof(double));  

  double *charge_old = (double *)calloc(n_grid, sizeof(double));

  double *potential = mPotentialProfile.GetPotential();

  double *pot_old = (double *)calloc(n_grid, sizeof(double));

  double phase;

  // Judge if outer self-consistent iteration is converged
  bool is_outer_converged = false;

  // Judge if inner poisson solver is converged
  bool is_inner_converged = true;

  // Count of the outer self-consistent iterations
  int count_outer_runs = 0;

  // Maximum potential change of the outer self-consistent iteration
  double max_pot_change;

  // Maximum charge change of the outer self-consistent iteration
  double max_charge_change; 

  // Records the smallest potential change of the outer self-consistent 
  // iteration (used for calling early stop)
  double min_max_pot_change = 9999.; 

  // Count of outer iterations whose potential change is large than previous
  // minimum value
  int count_early_stop = 0;


  // Print fixed charges
  if(mRank == 0)
  {
    cout << "Fixed charge profile" << endl;
    mFixedChargeProfile.Print();
  }


  // Runs outer simulation until converged
  while(is_outer_converged == false)
  {
    is_outer_converged = true;
    count_outer_runs++;

    // Stops when hitting maximum number of outer runs
    if(count_outer_runs > mSelfConsistentRunMax)
    {
      is_outer_converged = false;
      if(mRank == 0)
      {
        cout << "... Outer NEGF-Poisson solver not converged!" << endl;
      }
      break;
    }
    else
    {
      if(mRank == 0)
      {
        cout << endl << "Starting outer iteration number: " << count_outer_runs
          << " ..." << endl;
        cout << endl;
      }
    }

    // Records potential and charge from last iteration
    for(int i = 0; i < n_grid; ++i)
      pot_old[i] = potential[i];
    for(int i = 0; i < n_grid; ++i)
      charge_old[i] = charge[i];

    // Charge solver
    // Resets to be safe
    for(int i = 0; i < n_grid; ++i)
      charge[i] = 0.0;

    // Loops over assigned phases
    for(int i_phase = mIPhaseStart; i_phase < mIPhaseEnd; ++i_phase)
    {
      // Gets assigned phase value from index
      phase = mPhaseStep * double(i_phase);
      if(mRank == 0)
      {
        cout << "Computing phase point " << i_phase + 1 << " out of " <<
          mIPhaseEnd - mIPhaseStart << " ..." << endl;
        //        cout << "Current phase = " << phase / M_PI << " pi" << endl;
        cout << endl;
      }

      // Sets device Hamiltonians with current phase
      mChargeSolver.SetHamiltonians(phase);

      // Runs the charge solver to calculate charge for the given phase
      mChargeSolver.SolveCharge(mChargeProfile, mPotentialProfile);
    }

    // Collects charge from all processors by calling MPI_Allreduce()
    if(mRank == 0)
    {
      cout << "Waiting for other processors..." << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(mRank == 0)
    {
      cout << "All processors have finished!" << endl;
      cout << endl;
      cout << " Collecting data from all processors..." << endl;
      cout << endl;
    }
    MPI_Allreduce(charge, charge_sum, n_grid, MPI_DOUBLE, MPI_SUM,
        MPI_COMM_WORLD);

    // Averages over all phases to get the result
    for (int i = 0; i < n_grid; ++i)
    {
      charge_sum[i] /= (double)mNumKPoint;
      charge[i] = charge_sum[i];
    }

    // Updates charge profile with damping
    if(count_outer_runs > 1)
    {
      for(int i = 0; i < n_grid; ++i)
        charge[i] = charge_old[i] + (1.0 - mChargeDamping) * (charge[i]
            - charge_old[i]);
    }

    // Saves the updated charge profile
    mChargeProfile.SetCharge(charge);

    // Prints the updated charge profile
    if(mRank == 0)
    {
      cout << "Electron density profile updated from NEGF solver after damping\
        (A^-1):" << endl;
          mChargeProfile.Print();
    }
    if(mRank == 0)
    {
      cout << "Charge difference profile (A^-1):" << endl;
      for(int i = 0; i < n_grid; ++i)
      {
        cout << i + 1 << "\t" << std::setprecision(7) << charge_fixed[i] - charge[i]
          << endl;
      }
      cout << endl;
    }


    // Runs the poisson solver to get updated electrostatic potential
    mPoissonSolver.Run(mPotentialProfile, mChargeProfile, mFixedChargeProfile,
        is_inner_converged);

    // Stops self-consistent simulation if poisson solver is not converged
    if(is_inner_converged == false)
    {
      is_outer_converged = false;
      break;
    }

    // Updates potential profile with damping
    potential = mPotentialProfile.GetPotential();
    for(int i = 0; i < n_grid; ++i)
      potential[i] = pot_old[i] + (1.0 - mPotentialDamping) * (potential[i]
          - pot_old[i]);

    // Judges convergence of outer self-consistent iteration
    max_pot_change = 0.0;
    max_charge_change = 0.0;
    for(int i_pot = 0; i_pot < n_grid; ++i_pot)
    {
      if(max_pot_change < fabs(pot_old[i_pot] - potential[i_pot]))
        max_pot_change = fabs(pot_old[i_pot] - potential[i_pot]);
    }
    for(int i_pot = 0; i_pot < n_grid; ++i_pot)
    {
      if(max_charge_change < fabs(charge_old[i_pot] - charge[i_pot]))
        max_charge_change = fabs(charge_old[i_pot] - charge[i_pot]);
    }
    if(mRank == 0)
    {
      cout << "Maximum potential change = " << max_pot_change / (1.0
          - mPotentialDamping) << endl;
      cout << "Maximum charge distribution change = " << max_charge_change
        / (1.0 - mChargeDamping) << endl;
      cout << endl;
    }
    if(max_pot_change / (1.0 - mPotentialDamping) > mPotentialConv ||
        max_charge_change / (1.0 - mChargeDamping) > mPotentialConv * 1e1)
      is_outer_converged = false;

    // Check early stop
    if(max_pot_change < min_max_pot_change)
    {
      min_max_pot_change = max_pot_change;
      count_early_stop = 0;
    }
    else
      count_early_stop++;
    if(mRank == 0)
      cout << "Current early stop count is: " << count_early_stop << endl;      
    if(count_early_stop > mEarlyStop)
    {
      if(mRank == 0)      
        cout << "Early stop is reached!" << endl;
      break;
    }

    // Prints results if outer iteration is converged
    if(mRank == 0)
    {
      if(is_outer_converged == true)
      {
        cout << endl;
        cout << "Final result:" << endl;
        mPotentialProfile.Print();
        mChargeProfile.Print();
      }
    }
  }

  // Write results to file
  mOutput.WriteChargeToFile(mChargeProfile, mFixedChargeProfile,
      is_outer_converged);
  mOutput.WritePotentialToFile(mPotentialProfile, is_outer_converged);
}


  // Testing
  //int main()
  //{
  //  // Input parameters
  //  Parameters parameters;
  //  parameters.ParseInputFile();
  //  parameters.Print();  
  //
  //  Graphene* pGraphene = new Graphene();
  //  pGraphene->Init();
  //  
  //  MoS2* pMoS2 = new MoS2();
  //  pMoS2->Init();
  //
  //  SelfConsistentSolver central_region(parameters, graphene, mos2);
  //  central_region.Init();
  //
  //  return 0;
  //}
