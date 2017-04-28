# BootCMatch
Welcome to BootCMatch:
>     Bootstrap AMG based on Compatible Matching version 0.9
This is the software package described in the paper:
 P. D'Ambra, S. Filippone, P. Vassilevski
 BootCMatch: a software package for bootstrap AMG based on graph weighted matching
 submitted, 2017

It contains  the following components: 
1. A set of algorithms for building AMG hierarchies
2. A set of solvers including: 
  * Standalone multigrid
  * Flexible Conjugate Gradients preconditioned with AMG
3. A set of test programs.

The AMG hierarchies are built by using a graph matching approach. We support three different maximum weight matching algorithms:
1. The half-approximate Preis algorithm;
2. The approximate auction algorithm from [SPRAL](http://www.numerical.rl.ac.uk/spral/);
3. The exact matching from [HSL_MC64](http://www.hsl.rl.ac.uk/catalogue/hsl_mc64.html)

In the application of the AMG hierarchy (either as a solver or as a preconditioner) we provide a simple factorization, but we recommend to use our interface to [SuperLU](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/) for best performance.
To see the software in action we provide three test programs:
* testbootstrap     builds a set of AMG hierarchies
* testsolving       uses a set of AMG hierarchies as a linear system solver
* testkbootsolver   uses a set of AMG hierarchies as a preconditioner for the FCG method. 

At this time, you will have to edit the make.inc file by hand; when we get around to writing a configure script and a better  documentation we will declare version 1.0. 

For questions, please send an email to pasqua.dambra@cnr.it or salvatore.filippone@cranfield.ac.uk 

This project was partially supported by the NSF under award DMS-1619640.
This project was partially supported  by the EC under the Horizon 2020 Project Energy oriented Centre of Excellence for computing applications EoCoE, Project ID: 676629.
