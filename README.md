# BifurcationSnakesMilliKelvinSensitivity

This repository contains the C++ and Mathematica code used for the stochastic simulations and the parameter inference described in the paper "A bifurcation integrates information from many noisy ion channels and allows for milli-Kelvin thermal sensitivity in the snake pit organ" by Isabella Graf and Benjamin Machta.
For questions, please contact isabella.graf@yale.edu or benjamin.machta@yale.edu.



C++ code:

Thermosensing_Voltage_AP_Firing:
C++ code implementing the stochastic simulations of the time-discretized version of the full system, Eq. 2, in the paper.

Thermosensing_Voltage_AP_Firing_Simple_Model_Scaling:
C++ code implementing the stochastic simulations of the time-discretized version of the rescaled dynamics, Eq. 7, in the paper.

In both cases, parameters (and general settings) are defined in the main.cpp file, before compiling the code; comments in the code explain the meaning of the parameters and settings. 



Mathematica notebook:

compare_to_data_from_Bullock_Diecke.nb:
Mathematica notebook for the comparison to experimental data by Bullock and Diecke. The notebook implements the parameter inference and the corresponding best fits.
