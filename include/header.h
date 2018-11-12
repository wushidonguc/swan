// @author Wushi Dong

// header.h
// Purpose: Header file where we define global constants.

//#include <map>
//#include <string>
#include <complex>

#ifndef HEADER
#define HEADER

namespace constants
{
  const double BOLTZMANN_CONSTANT = 8.6173303e-5; /* eV * K^(-1) */
  const double ROOM_TEMPERATURE = 293.; /* K */
  const std::complex<double> I(0,1);
}

//namespace config
//{
//  // Map from MoS2 doping level to corresponding Fermi level
//  std::map<std::string,double> MoS2_doping_Fermi_level_map;
//  
//  MoS2_doping_Fermi_level_map['no doping'] = -2.1839;
//  MoS2_doping_Fermi_level_map["1e13"] = -2.1839 + 0.849121;
//  MoS2_doping_Fermi_level_map["5e13"] = -2.1839 + 0.972237;
//  MoS2_doping_Fermi_level_map["6e13"] = -2.1839 + 0.990981;
//  MoS2_doping_Fermi_level_map["7e13"] = -2.1839 + 1.0077;
//  MoS2_doping_Fermi_level_map["8e13"] = -2.1839 + 1.02242;
//  MoS2_doping_Fermi_level_map["1e14"] = -2.1839 + 1.04638;
//  MoS2_doping_Fermi_level_map["2e14"] = -2.1839 + 1.11189;
//  MoS2_doping_Fermi_level_map["3e14"] = -2.1839 + 1.15117;
//  MoS2_doping_Fermi_level_map["4e14"] = -2.1839 + 1.18351;
//  MoS2_doping_Fermi_level_map["5e14"] = -2.1839 + 1.21213;
//  MoS2_doping_Fermi_level_map["1e15"] = -2.1839 + 1.31826;
//  MoS2_doping_Fermi_level_map["1e13_low"] = -2.1839 + 0.879598;
//  MoS2_doping_Fermi_level_map["1e14_low"] = -2.1839 + 1.06333;
//  MoS2_doping_Fermi_level_map["2e14_low"] = -2.1839 + 1.12083;
//  MoS2_doping_Fermi_level_map["3e14_low"] = -2.1839 + 1.15578;
//}  

#endif
