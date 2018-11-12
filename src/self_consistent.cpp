// Swan: Self-consistent Wannier-function-based quantum transport solver
// Copyright (C) 2018 The University of Chicago
//
// This package is a free software.
// It is distributed under the xxx.
// The license text is found in the subfolder 'license' in the top folder.
// To request an official license document please write to the following address:
// 929 E 57th St, Chicago IL 60637, USA
//
// @author Wushi Dong
//
// self_consistent.cpp
//
// Purpose: Main function implementing the self-consistent simulation to obtain
// converged charge and electrostatic potential profile.

#include <iostream>
#include "mpi.h"

#include "Parameters.h"
#include "DeviceEdgeContact.h"
#include "SelfConsistentSolver.h"


int main(int argc, char * argv[])
{

  // Set up MPI environment
  int rank, nprocs;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	MPI_Comm cart_comm;

  // timers
  double startingTime = 0.;
  double finishingTime = 0.;


  // Print header
  if(rank == 0)
  {
    // Prints Logo
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "#####################################################" << std::endl;
    std::cout << "##                                                 ##" << std::endl;
    std::cout << "##                      Swan                       ##" << std::endl;
    std::cout << "##            Self-consistent simulaiton           ##" << std::endl;
    std::cout << "##                  by Wushi Dong                  ##" << std::endl;
    std::cout << "##                                                 ##" << std::endl;
    std::cout << "#####################################################" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    // Date and time
    time_t t = std::time(0);   // get time now
    tm* now = std::localtime(&t);
    std::cout << (now->tm_year + 1900) << '-' 
         << (now->tm_mon + 1) << '-'
         <<  now->tm_mday << ' '
         <<  now->tm_hour << ':'
         <<  now->tm_min << ':'
         <<  now->tm_sec
         << "\n";
    std::cout << std::endl;
  }

  
  // Print MPI environment information
  if(rank == 0)
    std::cout << "Number of running processors: " << nprocs << std::endl;


  // Get input parameters: create parameters with default values, parse input file, and print parameters
  if(rank == 0)
    std::cout << std::endl << "Parsing parameters ..." << std::endl << std::endl;
  Parameters parameters;
  parameters.ParseInputFile();
	if(rank == 0)
		parameters.Print();

   
  // Create materials
  if(rank == 0)
  {
    std::cout << "Creating materials ..." << std::endl << std::endl;
  }

  // Graphene
  if(rank == 0)
    std::cout << std::endl << "Graphene ..." << std::endl;
  const Graphene graphene(parameters);

  // MoS2
  if(rank == 0)
    std::cout << std::endl << "MoS2 ..." << std::endl;
  const MoS2 mos2(rank, parameters, graphene);
  

  // Set up the edge contact device
  if(rank == 0)
  {
    std::cout << std::endl << "Creating the edge contact device from materials ..." << std::endl << std::endl;
  }  
  const DeviceEdgeContact deviceEdgeContact(parameters, graphene, mos2);

  
  // Set up the self-consistent solver for the edge contact device
  if(rank == 0)
  {
    std::cout << std::endl << "Setting up the self-consistent solver ..." << std::endl << std::endl;
  }
  SelfConsistentSolver selfConsistentSolver(rank, nprocs, parameters, deviceEdgeContact);
  

  // Run the self-consistent transport simulation
  // Record starting time
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0)
	{
    startingTime = MPI_Wtime();
	}

  // Run
  if(rank == 0)
    std::cout << std::endl << "Running the self-consistent simulations ..." << std::endl << std::endl;
  selfConsistentSolver.Run();

  // Record finishing time
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0)
	{
    finishingTime = MPI_Wtime();
    std::cout << std::endl;
    std::cout << "Running time: " << ((int)(finishingTime - startingTime)) / 60 << " min " << ((int)(finishingTime - startingTime + 0.5)) % 60 << " sec" << std::endl << std::endl;
	}

  if(rank == 0)
    std::cout << "Completed!" << std::endl << std::endl;
  

  // Write results to file

  
  // Clean up

	MPI_Finalize();

  return 0;
}

