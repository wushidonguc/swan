// @author Wushi Dong

// TransportSolver.cpp

#include "TransportSolver.h"


TransportSolver::~TransportSolver(){}


void TransportSolver::Run()
{
  // Initialize
  // Number of grid point for the charge and potential profile
  const int n_grid = mConvergedPotentialProfile.GetNumGrid();
//  const int n_grid = 60;
//  std::cout << n_grid << endl;
  
  // Number of grid point for LDOS
  const int n_grid_LDOS = n_grid / 6 * 4 + 1;
  
  double * transmission = (double *)calloc(mEnergyRangeDense.mN,
      sizeof(double));
  double * transmission_sum = (double *)calloc(mEnergyRangeDense.mN,
      sizeof(double));
  
  double * LDOS = (double *)calloc(mEnergyRangeDense.mN * n_grid_LDOS,
      sizeof(double));
  double * LDOS_sum = (double *)calloc(mEnergyRangeDense.mN * n_grid_LDOS,
      sizeof(double));

  double phase;

  // Print converged potential
  if(mRank == 0)
  {
    cout << "The converged potential profile:" << endl;
    mConvergedPotentialProfile.Print();
  }

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
//    if(mRank == 0)
//      std::cout << "MARKER!" << endl;

    // Runs the charge solver to calculate charge for the given phase
    mChargeSolver.SolveTransmissionAndLDOS(transmission, LDOS, mConvergedPotentialProfile, mEnergyRangeDense);
//    if(mRank == 0)
//    std::cout << "MARKER3!!!" << endl;
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

  // Collects transmission
  MPI_Allreduce(transmission, transmission_sum, mEnergyRangeDense.mN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // Average over phases
  for (int i = 0; i < mEnergyRangeDense.mN; ++i)
  {
    transmission_sum[i] /= (double)mNumKPoint;
  }

  // Collects LDOS
  MPI_Allreduce(LDOS, LDOS_sum, mEnergyRangeDense.mN * n_grid_LDOS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // Average over phases
  for (int i = 0; i < mEnergyRangeDense.mN * n_grid_LDOS; ++i)
  {
    LDOS_sum[i] /= (double)mNumKPoint;
  }

  // Write results to file
  mOutput.WriteTransmissionToFile(transmission_sum, mEnergyRangeDense);
  mOutput.WriteLDOSToFile(LDOS_sum, n_grid_LDOS, mEnergyRangeDense);
}

