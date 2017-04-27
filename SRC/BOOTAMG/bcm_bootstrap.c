/* 
                BootCMatch 
     Bootstrap AMG based on Compatible Matching version 0.9
    (C) Copyright 2017
                       Pasqua D'Ambra         IAC-CNR, IT
                       Salvatore Filippone    Cranfield University, UK
                       Panayot S. Vassilevski Portland State University, OR USA
 
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
    1. Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above copyright
       notice, this list of conditions, and the following disclaimer in the
       documentation and/or other materials provided with the distribution.
    3. The name of the BootCMatch group or the names of its contributors may
       not be used to endorse or promote products derived from this
       software without specific written permission.
 
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE BootCMatch GROUP OR ITS CONTRIBUTORS
  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.
 
*/

#include "bcm_bootamg.h"

/* ************************************************************************
*
*  function for generating multiple AMG hierarchy to build a composite solver 
*  having a desired convergence rate
*
************************************************************************** */

bcm_BootAMG *
bcm_Bootstrap(bcm_BootAMGBuildData *bootamg_data, bcm_AMGApplyData *amg_cycle)

{

  double desired_ratio;
  int max_hrc;
  int ierr;

  max_hrc=bcm_BootAMGBuildDataMaxHrc(bootamg_data);
  desired_ratio=bcm_BootAMGBuildDataDesRatio(bootamg_data);

  bcm_BootAMG *boot_amg;
  boot_amg=bcm_BootAMGCreate(max_hrc);
  bcm_BootAMGInitialize(boot_amg);
  bcm_AMGHierarchy **Harray;
  Harray=bcm_BootAMGHarray(boot_amg); 

  double conv_ratio=bcm_BootAMGEstRatio(boot_amg);

  int num_hrc=0, i;

  bcm_AMGBuildData *amg_data;
  amg_data=bcm_BootAMGBuildDataCompData(bootamg_data);

  while(conv_ratio > desired_ratio & num_hrc < max_hrc)
    {
      bcm_AMGHierarchyDestroy(Harray[num_hrc]);	 
      Harray[num_hrc]=bcm_AdaptiveCoarsening(amg_data);
      num_hrc++;
      bcm_BootAMGNHrc(boot_amg)=num_hrc;
      printf("Built new hierarchy. Current number of hierarchies:%d %d\n",num_hrc,max_hrc);
      if(max_hrc > 1)
	{
	  ierr=bcm_InnerIterations(bootamg_data, boot_amg, amg_cycle);
	  conv_ratio=bcm_BootAMGEstRatio(boot_amg);
	  printf("Convergence ratio of the current composite solver =%e\n",conv_ratio);
	}
    }
  ierr=bcm_BootAMGHarrayDestroy(max_hrc-num_hrc, num_hrc, Harray);

  return boot_amg;
}


/* ************************************************************************
*
*  function for applying a composite AMG for a fixed number of iterations
*  to homogeneous system associated to the given matrix, needed to estimate 
*  solver convergence ratio
*
************************************************************************** */

int
bcm_InnerIterations(bcm_BootAMGBuildData *bootamg_data, bcm_BootAMG *boot_amg, bcm_AMGApplyData *amg_cycle)

{
  int ierr=0;
  double conv_ratio, normold, normnew, alpha;
  int iter, n_hrc, solver_it, solver_type;
  int i,k;
  int num_lev;

  bcm_Vector **Xtent, **RHS;

  conv_ratio=bcm_BootAMGEstRatio(boot_amg);
   
  bcm_AMGHierarchy **Harray;
  Harray=bcm_BootAMGHarray(boot_amg);

  n_hrc=bcm_BootAMGNHrc(boot_amg);

  /* data related to single AMG */
  bcm_AMGBuildData *amg_data;
  amg_data=bcm_BootAMGBuildDataCompData(bootamg_data);

  /* data problem (original matrix and current smooth vector) related to single AMG */
  bcm_CSRMatrix *A;
  bcm_Vector *w;
  A=bcm_AMGBuildDataCSRMatrix(amg_data);
  w=bcm_AMGBuildDataSmoothVector(amg_data);
  int nsize=bcm_VectorSize(w);

  solver_it=bcm_BootAMGBuildDataCompIt(bootamg_data);
  solver_type=bcm_BootAMGBuildDataCompType(bootamg_data);

  bcm_Vector *x; /* current solution: at the end it will point to the new smooth vector */
  bcm_Vector *rhs; /* rhs */

  x=bcm_VectorCreate(nsize);
  rhs=bcm_VectorCreate(nsize);
  bcm_VectorInitialize(x);
  bcm_VectorInitialize(rhs);
  bcm_VectorSetConstantValues(rhs,0.0);
  bcm_VectorSetRandomValues(x,1.0);

  normold=bcm_VectorANorm(A,x);

  bcm_CSRMatrix **A_array;

  /* start iterative cycle */
  for(iter=1; iter<=solver_it; iter++)
    {
      ierr = bcm_PrecApply(bootamg_data, boot_amg, amg_cycle, rhs, x);
      normnew=bcm_VectorANorm(A,x);
      conv_ratio=normnew/normold;
      normold=normnew;
    }

  /* update smooth vector */
  alpha=1.0/normnew;

  printf("current smooth vector A-norm=%e\n",normnew);

  bcm_VectorScale(alpha,x);
  bcm_VectorCopy(x,w); /* update the smooth vector for the possible new bootstrap step */

  bcm_BootAMGEstRatio(boot_amg)=conv_ratio;

  bcm_VectorDestroy(x);
  bcm_VectorDestroy(rhs);

  return ierr;
}
