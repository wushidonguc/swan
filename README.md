# Swan
An open-source C++ software for nanoscale quantum electron transport simulations 

**Swan** (**S**elf-consistent **wan**nier-function-based quantum transport solver) is an open-source C++ software suitable for large-scale atomistic simulations of electronic structure and transport properties in nano-devices.
By using a Wannier function basis (as implemented in the [Wannier90](http://wannier.org) package) to accurately describe the electronic bands, our code is able to efficiently model device structures with first-principles accuracy at a minimal cost of tight-binding calculations.
It couples the Keldysh non-equilibrium Green's function formalism and Poisson solver to generate the inhomogeneous charge densities self-consistently with the electrostatic potential profile for the simulated device region. The parallel implementation of the code uses the standard Message Passing Interface (MPI).

For more details, please refer to the related paper "Quantum electron transport in ohmic edge contacts between two-dimensional materials" ([https://arxiv.org/abs/1811.02135](https://arxiv.org/abs/1811.02135)).

Author: Wushi Dong (dongws@uchicago.edu) of the Physics Department at The University of Chicago. (Advisor: Peter B. Littlewood)

The author would like to acknowledge the C++ linear algebra library Armadillo.

# Installation

Please install the Armadillo package via the download page:

[http://arma.sourceforge.net/download.html](http://arma.sourceforge.net/download.html)

Then install Swan by running the installation file provided in the package:

```
sh install.sh
```

# Example run

To demonstrate the usage of our code, we provide an example run which simulates the 2D edge contact structure studied in our paper ([https://arxiv.org/abs/1811.02135](https://arxiv.org/abs/1811.02135)). To run this example:

```
cd example_run/
sh run.sh
```
