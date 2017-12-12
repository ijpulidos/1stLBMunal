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
