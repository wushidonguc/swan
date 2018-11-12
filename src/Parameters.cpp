// @author Wushi Dong

// Parameters.cpp

#include "Parameters.h"

/* Parameters class */

Parameters::Parameters()
{
  mTemperature = constants::ROOM_TEMPERATURE;
  mNumPrincipalLayerGraphene = 6;
  mNumPrincipalLayerMoS2 = 4;
  mFermiLevelGr = 0.0;
  mDoping = "1e14";
  mDistanceGrS2 = 1.5;
  mHoppingCS = -1.0;
  mEnergyMin = -9.0;
  mEnergyMax = 5.0;
  mEnergyStep = 0.01;
//  mEnergyMinDense = -9.0;
//  mEnergyMaxDense = 5.0;
//  mEnergyStepDense = 0.002;
  mNumKPoint = 32;
  mVoltageBias = 0.0;
  mSurfaceGreensFunctionConv = 1.0e-9;
  mPotentialConv = 5.0e-4;
  mSelfConsistentRunMax = 2000;
  mPoissonRunMax = 200;
  mPotentialDamping = 0.1;
  mChargeDamping = 0.1;
  mPotentialDampingPoisson = 0.9;
  mEarlyStop = 30;
  mPotentialInitializationMethod = "linear";
  mHoppingFile = "mos2_hr.txt";
  mFixedChargeFileGr = "fixed_charge_Gr.txt";
  mFixedChargeFileMoS2 = "fixed_charge_MoS2_1e14.txt";
}
Parameters::~Parameters(){}

void Parameters::ParseInputFile()
{
  std::ifstream infile;
  infile.open("parameters.txt");

  std::string name, value;
  while( infile >> name >> value )
  {
    if(name[0] == '#') continue;
    if(name[0] == '\n') continue;

    else if(name == "temperature")
      mTemperature = stod(value);

    else if(name == "n_Gr")
      mNumPrincipalLayerGraphene = stoi(value);

    else if(name == "n_MoS2")
      mNumPrincipalLayerMoS2 = stoi(value);

    else if(name == "E_Fermi_Gr")
      mFermiLevelGr = stod(value);    

    else if(name == "doping")
      mDoping = value;

    else if(name == "distance_Gr_S2")
      mDistanceGrS2 = stod(value);    

    else if(name == "hopping_CS")
      mHoppingCS = stod(value);

    else if(name == "energy_min")
      mEnergyMin = stod(value);    

    else if(name == "energy_max")
      mEnergyMax = stod(value);    

    else if(name == "energy_step")
      mEnergyStep = stod(value);
    
    else if(name == "n_k_point")
      mNumKPoint = stoi(value);

    else if(name == "voltage_bias")
      mVoltageBias = stod(value);
    
    else if(name == "surface_conv")
      mSurfaceGreensFunctionConv = stod(value);
    
    else if(name == "pot_conv")
      mPotentialConv = stod(value);
    
    else if(name == "sc_run_max")
      mSelfConsistentRunMax = stoi(value);

    else if(name == "poisson_run_max")
      mPoissonRunMax = stoi(value);
    
    else if(name == "pot_damping")
      mPotentialDamping = stod(value);
    
    else if(name == "charge_damping")
      mChargeDamping = stod(value);

    else if(name == "pot_damping_poisson")
      mPotentialDampingPoisson = stod(value);

    else if(name == "early_stop")
      mEarlyStop = stoi(value);

    else if(name == "pot_init")
      mPotentialInitializationMethod = value;

    else if(name == "hopping_file")
      mHoppingFile = value;

    else if(name == "fixed_charge_Gr")
      mFixedChargeFileGr = value;

    else if(name == "fixed_charge_MoS2")
      mFixedChargeFileMoS2 = value;
  }
  
  infile.close();
  
}

// Print parameters
void Parameters::Print() const
{
  std::cout << "**************************PARAMETERS*************************" << std::endl;
  std::cout << "temperature = " << mTemperature << " K" << std::endl;
  std::cout << "n_Gr = " << mNumPrincipalLayerGraphene << std::endl;
  std::cout << "n_MoS2 = " << mNumPrincipalLayerMoS2 << std::endl;
  std::cout << "E_Fermi_Gr = " << mFermiLevelGr << std::endl;
  std::cout << "doping = " << mDoping << std::endl;
  std::cout << "distance_Gr_S2 = " << mDistanceGrS2 << " Ang" << std::endl;
  std::cout << "hopping_CS = " << mHoppingCS << " eV" << std::endl;
  std::cout << "energy_min = " << mEnergyMin << " eV" << std::endl;
  std::cout << "energy_max = " << mEnergyMax << " eV" << std::endl;
  std::cout << "energy_step = " << mEnergyStep << " eV" << std::endl;
  std::cout << "n_k_point = " << mNumKPoint << std::endl;
  std::cout << "voltage_bias = " << mVoltageBias << " V" << std::endl;
  std::cout << "surface_conv = " << mSurfaceGreensFunctionConv << " eV" << std::endl;
  std::cout << "pot_conv = " << mPotentialConv << " V" << std::endl;
  std::cout << "sc_run_max = " << mSelfConsistentRunMax << std::endl;
  std::cout << "poisson_run_max = " << mPoissonRunMax << std::endl;
  std::cout << "pot_damping = " << mPotentialDamping << std::endl;
  std::cout << "charge_damping = " << mChargeDamping << std::endl;
  std::cout << "pot_damping_poisson = " << mPotentialDampingPoisson << std::endl;
  std::cout << "early_stop = " << mEarlyStop << std::endl;
  std::cout << "pot_init = " << mPotentialInitializationMethod << std::endl;
  std::cout << "hopping_file = " << mHoppingFile << std::endl;
  std::cout << "fixed_charge_Gr = " << mFixedChargeFileGr << std::endl;
  std::cout << "fixed_charge_MoS2 = " << mFixedChargeFileMoS2 << std::endl;
  std::cout << "*************************************************************" << std::endl;
}


//// Testing
//int main()
//{
//  Parameters parameters;
//  parameters.ParseInputFile();
//  parameters.Print();
// 
//  return 0;
//}
