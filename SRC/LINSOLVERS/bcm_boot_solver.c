/*
                BootCMatch
     Bootstrap AMG based on Compatible weighted Matching, version 0.9
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

#include "bcm_linsolvers.h"

/* ************************************************************************
*
*  function for applying a composite AMG, built by the bootstrap process,
*  as a solver on its own, until convergence
*
************************************************************************** */

int bcm_boot_solver(bcm_BootAMGBuildData *bootamg_data, bcm_BootAMG *boot_amg,
		bcm_AMGApplyData *amg_cycle,
		bcm_Vector *rhs, bcm_Vector *x,  int solver_it, double RTOL)

{
  int ierr=0;
  double normold, normnew, alpha;
  int iter, n_hrc, solver_type;
  int i,k;
  int num_lev;

  bcm_Vector **Xtent, **RHS;

  bcm_AMGHierarchy **Harray;
  Harray=bcm_BootAMGHarray(boot_amg);
  n_hrc=bcm_BootAMGNHrc(boot_amg);

  /* data related to single AMG */
  bcm_AMGBuildData *amg_data;
  amg_data=bcm_BootAMGBuildDataCompData(bootamg_data);

  /* starting data */
  bcm_CSRMatrix *A;
  A=bcm_AMGBuildDataCSRMatrix(amg_data);
  int nsize=bcm_CSRMatrixNumRows(A);
  solver_type=bcm_BootAMGBuildDataCompType(bootamg_data);

  /* compute initial residual */

  bcm_Vector *res;
  res=bcm_VectorCreate(nsize);
  bcm_VectorInitialize(res);
  bcm_VectorCopy(rhs,res);
  bcm_CSRMatrixMatvec(-1.0,A,x,1.0,res);

  normold=bcm_VectorANorm(A,res);

  /* start iterative cycle */
  double resinit=normold;
  double relres=1.0;
  int j=1;

  bcm_CSRMatrix **A_array;

  while(relres >RTOL & j<=solver_it)
    {
      ierr = bcm_PrecApply(bootamg_data, boot_amg, amg_cycle, rhs, x);

      bcm_VectorCopy(rhs,res);
      bcm_CSRMatrixMatvec(-1.0,A,x,1.0,res);
      normnew=bcm_VectorANorm(A,res);
      relres=normnew/resinit;
      normold=normnew;
      printf("relative residual at iteration it=%d, %e\n",j,relres);
      j++;
    }
  bcm_VectorDestroy(res);
  if(j> solver_it) printf("Warning: Accuracy Requests was not satisfied!\n");
  printf("Number of Iterations =%d \n", j-1);
  printf("current residual A-norm =%e \n", normnew);
  printf("relative residual in A-norm=%e \n", relres);

  return ierr;
}
