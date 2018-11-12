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
// transport.cpp
//
// Purpose: Main function computing the transmission spectrum and Local Density
// of States (LDOS) given the converged electrostatic potential profile from the
// self-consistent simualtion.

#include <iostream>
#include "mpi.h"

#include "Parameters.h"
#include "DeviceEdgeContact.h"
#include "TransportSolver.h"


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
    std::cout << endl;
    std::cout << endl;
    std::cout << endl;
    std::cout << "#####################################################" << endl;
    std::cout << "##                                                 ##" << endl;
    std::cout << "##                      Swan                       ##" << endl;
    std::cout << "##        (Calculate transmission and LDOS)        ##" << endl;
    std::cout << "##                  by Wushi Dong                  ##" << endl;
    std::cout << "##                                                 ##" << endl;
    std::cout << "#####################################################" << endl;
    std::cout << endl;
    std::cout << endl;

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
    std::cout << endl;
  }

  
  // Print MPI environment information
  if(rank == 0)
    std::cout << "Number of running processors: " << nprocs << endl;


  // Get input parameters: create parameters with default values, parse input file, and print parameters
  if(rank == 0)
    std::cout << endl << "Parsing parameters ..." << endl << endl;
  Parameters parameters;
  parameters.ParseInputFile();
	if(rank == 0)
		parameters.Print();

   
  // Create materials
  if(rank == 0)
  {
    std::cout << endl << "Creating materials ..." << endl;
  }

  // Graphene
  if(rank == 0)
  {
    std::cout << "Graphene ..." << endl;
  }
  const Graphene graphene(parameters);

  // MoS2
  if(rank == 0)
  {
    std::cout << "MoS2 ..." << endl;
  }
  const MoS2 mos2(rank, parameters, graphene);
  

  // Set up the edge contact device
  if(rank == 0)
  {
    std::cout << endl << "Using created materials to make the edge contact device  ..." << endl;
  }  
  const DeviceEdgeContact deviceEdgeContact(parameters, graphene, mos2);

  
  // Set up the transport solver for the edge contact device
  if(rank == 0)
  {
    std::cout << endl << "Setting up the transport solver ..." << endl;
  }
  TransportSolver transportSolver(rank, nprocs, parameters, deviceEdgeContact);
  

  // Run the transport calculation
  // Record starting time
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0)
	{
    startingTime = MPI_Wtime();
	}

  // Run
  if(rank == 0)
    std::cout << endl << "Running the transport calculations ..." << endl << endl;
  transportSolver.Run();

  // Record finishing time
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0)
	{
    finishingTime = MPI_Wtime();
    std::cout << endl;
    std::cout << "Running time: " << ((int)(finishingTime - startingTime)) / 60 << " min " << ((int)(finishingTime - startingTime + 0.5)) % 60 << " sec" << endl << endl;
	}

  if(rank == 0)
    std::cout << "Completed!" << endl << endl;
  
 
  // Clean up

	MPI_Finalize();

  return 0;
}

