// @author Wushi Dong

// Output.cpp

#include "Output.h"


Output::~Output(){}


void Output::WriteChargeToFile(const ChargeProfile charge_profile, const
    FixedChargeProfile fixed_charge_profile, const bool &rIsOuterConverged)
const
{
  if(mRank == 0)
  {
    double *charge = charge_profile.GetCharge();
    double *charge_fixed = fixed_charge_profile.GetCharge();
    int n_grid = charge_profile.GetNumGrid();
    std::ofstream myfile;

    // If outer run is converged
    if(rIsOuterConverged == true)
    {
      myfile.open ("out/charge_conv.txt");
      for (int i = 0; i < n_grid; ++i)
      {
        myfile << i + 1 << "\t" << std::setprecision(7) << charge_fixed[i]
          - charge[i] << "\t" << std::setprecision(7) << -charge[i] << endl;
      }
      myfile.close();
    }
    // If outer run is not converged
    else
    {
      std::ofstream myfile;
      myfile.open ("out/charge_difference_not_conv.txt");
//      myfile << "# max_charge_change = " << max_charge_change << endl;
      for (int i = 0; i < n_grid; ++i)
      {
        myfile << i + 1 << "\t" << std::setprecision(7) << charge_fixed[i]
          - charge[i] << "\t" << std::setprecision(7) << -charge[i] << endl;

      }
      myfile.close();
    }
  }
}


void Output::WritePotentialToFile(const PotentialProfile
    potential_profile, const bool &rIsOuterConverged) const
{
  if(mRank == 0)
  {
    double *potential = potential_profile.GetPotential();
    int n_grid = potential_profile.GetNumGrid();
    std::ofstream myfile;

    // If outer run is converged
    if(rIsOuterConverged == true)
    {
      myfile.open ("out/potential_conv.txt");
      for (int i = 0; i < n_grid; ++i)
        myfile << i + 1 << "\t" << std::setprecision(7) << potential[i] << endl;
      myfile.close();
    }
    // If outer run is not converged
    else
    {
      myfile.open ("out/potential_not_conv.txt");
//      myfile << "# max_pot_change = " << max_pot_change << endl;
//      myfile << "# early_stop is " << std::boolalpha << bool(early_stop_count > early_stop) << endl;
      for (int i = 0; i < n_grid; ++i)
        myfile << i + 1 << "\t" << std::setprecision(7) << potential[i] << endl;
      myfile.close();
    }
  }
}


void Output::WriteTransmissionToFile(double * transmission, EnergyRange energyRangeDense) const
{
  if(mRank == 0)
  {
    std::ofstream myfile;
    myfile.open ("out/transmission.txt");
  
    // If outer run is converged
    for (int i = 0; i < energyRangeDense.mN; ++i)
    {
      myfile << energyRangeDense.mStart + (double(i) + 0.5) * energyRangeDense.mSpacing << "\t" <<
        std::setprecision(7) << transmission[i] << endl;
    }
      myfile.close();
  }
}

// TODO(dongws@uchicago.edu): print positions of LDOS grids
void Output::WriteLDOSToFile(double * LDOS, int n_grid_LDOS, EnergyRange energyRangeDense) const
{
  if(mRank == 0)
  {
    std::ofstream myfile;
    myfile.open ("out/LDOS.txt");
    for(int i = 0; i < n_grid_LDOS; ++i)
    {
      for (int j = 0; j < energyRangeDense.mN; ++j)
      {
        myfile << i << "\t" << energyRangeDense.mStart + (double(j) + 0.5) * energyRangeDense.mSpacing << "\t" << std::setprecision(7) << LDOS[i * energyRangeDense.mN + j] << endl;
      }
      myfile << endl;
    }
    myfile.close();
  }
}

