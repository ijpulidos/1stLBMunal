# 1st course on Lattice-Boltzmann methods - National University of Colombia - Dic2017

This repository will save notes from the 1st course on Lattice-Boltzmann methods at National 
University of Colombia.

## Using this repository

To use this you just have to clone it using the following command

    git clone https://github.com/ijpulidos/1stLBMunal.git

and then change to the new created directory

    cd 1stLBMunal
    
where you can find the material distributed in different topics and dates.

## LB method for waves - 2017-12-11

This first session was about creating a Lattice-Boltzmann method implementation for simple 
scalar waves. Go to directory `LBwaves-2017-12-11` and just compile the code with

    g++ LB_waves.cpp
    
which creates an executable file `a.out` and then run it using

    ./a.out
    
This will create a file called `ondas.dat` (sorry for the spanglish ;) with the data to be 
plotted using gnuplot, by running the `gnuplot` command and inside the gnuplot terminal, the
following
    
    set pm3d
    unset surface
    splot "ondas.dat"

## LB method for Poisson eq. - 2017-12-12

The idea in this session was completing the code for the LBM for the Poisson equation.

Completing/implementing density charges and initial conditions and finding a power law for
the error coefficient epsilon respect to the size of the system.

**Solution**

The solution to the exercise is in both `LBPoisson_solution.ccp` and my solution in 
`LBPoisson_exercise.cpp`. A plotting routine to plot the RMS error in terms of the size of
the lattice is in the python script `plotting.py`

## LBM for fluid dynamics - 2017-12-13

LBM for basic fluid dynamics, relatively low Reynolds numbers.

This would create a file called "fluid.dat" with x, y, vx, vy data in it. You
can plot the data using `gnuplot` with the following commands for heat map:

  set pm3d
  unset surface
  set view map
  splot "fluid.dat"

or to see arrows

  plot "fluid.dat" with vector
    
