## Political Economy of Taxation in an Overlapping-Generations Economy
*Review of Economic Dynamics*, 2008, vol. 11, n.1, pp.18-43.

This code is written for Matlab 6.

Main file: [olgsolve.m](code/olgsolve.m)

Subroutine needed for the computation of the threat point, under endogenous and exogenous prices: [kstexo.m](code/kstexo.m), [kstendo.m](code/kstendo.m)

Objective function in solving for the political-economic equilibrium: [nashfunexo.m](code/nashfunexo.m), [nashfunendo.m](code/nashfunendo.m)

Contraints in solving for the political-economic equilibrium: [nashconexo.m](code/nashconexo.m), [nashconendo.m](code/nashconendo.m)

Subroutine needed for the computation of the steady state: [steady.m](code/steady.m)

Main file in the case of no government spending: [olgsolvenog.m](code/olgsolvenog.m)

Objective function in solving for the political-economic equilibrium with no government spending: [nashfunendonog.m](code/nashfunendonog.m)

Constraints in solving for the political-economic equilibrium with no government spending: [nashconendonog.m](code/nashconendonog.m)

Note: This code also requires the subroutine for the computation of the threat point, which is the same as above ([kstendo.m](code/kstendo.m)). 