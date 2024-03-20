## Politics and Efficiency of Separating Capital and Ordinary Government Budgets (with Thomas J. Sargent)
*Quarterly Journal of Economics*, 2006, vol. 121 n. 4,  pp. 1167-1210.

### Code for the calibration of section 4

This code is written for Matlab 6.

Main file that loads demographics and depreciation for all the different cases: [multi.m](code/demographics/multi.m)

Script file that uses the parameter values from multi.m to generate all the results: [hetprob.m](code/demographics/hetprob.m)

Files that contain demographic data needed as input:

Mortality by age in the U.S. 2000: [usdeathsnow.txt](code/demographics/usdeathsnow.txt)

Age structure of the U.S. population in 2000: [usagenow.txt](code/demographics/usagenow.txt)

Mortality by age in the U.S. 1880: [usdeaths1880.txt](code/demographics/usdeaths1880.txt)

Age structure of the U.S. population in 1880: [usage1880.txt](code/demographics/usage1880.txt)

Mortality by age in the U.S. 1880 (males only): [usdm1880.txt](code/demographics/usdm1880.txt)

Age structure of the U.S. male population in 1880: [usam1880.txt](code/demographics/usam1880.txt)

Population growth by state between 1990 and 2000: [popgrowthallstates.txt](code/demographics/popgrowthallstates.txt)

Table with migration data by state and age in 2000 (after an appropriate reshaping, described in multi.m, the 6th column contains emigration between 1995 and 2000, the 1st column contains the number of people present in 2000, the 7th column contains net migration between 1995 and 2000, and the 8th column contains immigration from abroad between 1995 and 2000): [migrationallstates.txt](code/demographics/migrationallstates.txt)

Age structure of the population by state in 1990: [pop1990allstates.txt](code/demographics/pop1990allstates.txt)

Age structure of the population by state in 2000: [pop2000allstates.txt](code/demographics/pop2000allstates.txt)

### Code for appendix 3

The code for the log case is in Fortran. The Fortran code comprises 2 files that vary by experiment, and 5 that are common to all experiments. 

Experiment 1: [modelsetup.f90](code/0flatillinois/modelsetup.f90), [br18bclsj.f90](code/0flatillinois/br18bclsj.f90). The first file provides parameter values, the second file is the main file that calls all other subroutines.

Experiment 2: [modelsetup.f90](code/1flatmaxconflict/modelsetup.f90), [br18bclsj.f90](code/1flatmaxconflict/br18bclsj.f90)

Experiment 2 with theta_0=0.4: [modelsetup.f90](code/2flatlessconflict/modelsetup.f90), [br18bclsj.f90](code/2flatlessconflict/br18bclsj.f90)

Experiment 2 with theta_0=0.4 and an increasing consumption profile:  [modelsetup.f90](code/3incrlessconflict/modelsetup.f90), [br18bclsj.f90](code/3incrlessconflict/br18bclsj.f90)

Common files:

[mainloop.f90](code/common18/mainloop.f90) contains the main subroutines that iterate to convergence of the political-economic equilibrium and compute the steady state

[chebyshev.f90](code/common18/chebyshev.f90) contains many subroutines used for Chebyshev interpolation

[fcn.f90](code/common18/fcn.f90) function providing the residual from the first-order conditions that must hold in a (differentiable, Markov) political-economic equilibrium 

[lsjac.f90](code/common18/lsjac.f90) contains the Jakobian of the function fcn

[writelog.f90](code/common18/writelog.f90) subroutines to print and visualize the output

The code for the linear case is in Matlab. It is composed of a singlle file for each of the experiments; the only difference across experiments are the parameter values.

Experiment 1: [hetprob.m](code/0flatillinois/hetprob.m)

Experiment 2: [hetprob.m](code/1flatmaxconflict/hetprob.m)

Experiment 2 with theta_0=0.4 (it applies to both constant and increasing consumption profile, see paper): [hetprob.m](code/2flatlessconflict/hetprob.m)