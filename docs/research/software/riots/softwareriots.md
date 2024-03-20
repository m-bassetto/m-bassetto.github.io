## Tax Riots (joint with Christopher Phelan)

*Review of Economic Studies*, 2008, vol. 75, n.3, pp.649-669.

The code for this project is written for Mathematica.

- Main code, that reproduces the graphs of the numerical example in the paper: [eqnsolver052807.nb](code/eqnsolver052807.nb)
- Code for the degenerate case. While the previous file works fine even for a degenerate case, this file contains computations for the case in which the optimal ex ante plan is perturbed on a small-probability event (to check the claim that, in an economy with 1,000,000 people changing what happens on a 1% probability event is well more than needed to obtain that truth telling is strictly interim dominant): [degenerate2.nb](code/degenerate2.nb)
- Code that includes a computation of the maximin (it is slightly less general than the main code, since it assumes full support; by default, it is set to use a truncated normal distribution for f, but this can be changed): [eqnsolvernormal2.nb](code/eqnsolvernormal2.nb)
- Code for computing mixed-strategy equilibria when resources seized from tax evaders can be redistributed: [redistributing.nb](code/redistributing.nb)
- Code for the case in which taxes on households that report low income are constrained to be the same (assuming that they are not discovered to be cheaters): [taulowthesame.nb](code/taulowthesame.nb)
- Code for a CRRA Mirrlees economy with a degenerate distribution f: [mirrleesbothworkdegenerate.nb](code/mirrleesbothworkdegenerate.nb)
- Code for a CRRA Mirrlees economy with uniform distribution f: [mirrleesbothwork.nb](code/mirrleesbothwork.nb)
- Code for a CARA Mirrlees economy with degenerate distribution f: [mirrleescarabothworkdegenerate.nb](code/mirrleescarabothworkdegenerate.nb)
- Code for a Green economy with a degenerate distribution f: [green.nb](code/green.nb)
