# BootCMatch
Welcome to BootCMatch:
>     Bootstrap AMG based on Compatible weighted Matching, version 1.0.

Pasqua D'Ambra, Institute for Applied Computing (IAC) of the Italian National Research Council (CNR), Italy

Panayot S. Vassilevski, Portland State University, OR USA

This is the software package described in the paper:

P. D'Ambra, S. Filippone, P. S. Vassilevski: *BootCMatch: a software package for bootstrap AMG based on graph weighted matching,*
ACM Transactions on Mathematical Software, Vol. 44, 39:1-39:25, 2018

It contains the following components: 
1. A set of algorithms for building AMG hierarchies
2. A set of solvers including: 
  * Standalone multigrid
  * Flexible Conjugate Gradients preconditioned with AMG
3. A set of test programs.

The AMG hierarchies are built by using a weighted graph matching approach. We support four different maximum weight matching algorithms:
1. The half-approximate algorithm (Preis, 1999);
2. The approximate auction algorithm from [SPRAL] (http://www.numerical.rl.ac.uk/spral/);
3. The optimal matching from [HSL_MC64] (http://www.hsl.rl.ac.uk/catalogue/hsl_mc64.html);
4. The half-approximate Suitor algorithm (Manne, Halappanavar, 2014); 

In the application of the AMG hierarchy (either as a solver or as a preconditioner) we provide a simple factorization, but we recommend to use our interface to [SuperLU] (http://crd-legacy.lbl.gov/~xiaoye/SuperLU/) for best performance.
To see the software in action we provide three test programs:
* testbootstrap     builds a set of AMG hierarchies
* testsolving       uses a set of AMG hierarchies as a linear system solver
* testkbootsolver   uses a set of AMG hierarchies as a preconditioner for the FCG method. 

At this time, you will have to edit the make.inc file by hand to suit your computing environment. 

Contributors for version 0.9:
       Salvatore Filippone. 

For questions, please send an email to pasqua.dambra@cnr.it or open an
issue on our Git page: https://github.com/bootcmatch/BootCMatch

This project was partially supported by the NSF under award DMS-1619640.
This project was partially supported  by the EC under the Horizon 2020
Project: *Energy oriented Centre of Excellence for computing
applications EoCoE*, Project ID: 676629. 
