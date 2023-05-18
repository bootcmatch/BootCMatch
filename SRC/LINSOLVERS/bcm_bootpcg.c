/*
                BootCMatch
     Bootstrap AMG based on Compatible weighted Matching, version 1.0
    (C) Copyright 2017
                       Pasqua D'Ambra         IAC-CNR, IT
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

/* This function applies a Preconditioned Flexible CG (FCG) Krylov method (only
 * 1 direction is used at each iteration), where the preconditioner is the 
 * first hierarchy built by the bootstrap process. 
 * N.B. Only 1 AMG hierarchy based on compatible weighted matching is applied  */

#include "bcm_linsolvers.h"

int 
bcm_bootpcg(bcm_Vector *x_vector, bcm_Vector *rhs_vector, bcm_BootAMGBuildData *bootamg_data,
	    bcm_BootAMG * boot_amg, bcm_AMGApplyData *amg_cycle, int precon, int max_iter,
	    double rtol, int *num_iter, double *timetot)

{

  int ierr=0, i;
  double rhs_norm ,delta0, delta_old, delta, eps = DBL_EPSILON; 
  bcm_Vector **d_vector, *v_vector, *w_vector; 
  int num_dofs;

  int iter=0, idx; 
  double tau, tau1, alpha, beta;

  double l2_norm;
  double   zero=0.0, one=1.0, minusone=-1.0;

  /* matrix data */
  bcm_AMGBuildData *amg_data;
  amg_data=bcm_BootAMGBuildDataCompData(bootamg_data);

  bcm_CSRMatrix *A;
  A=bcm_AMGBuildDataCSRMatrix(amg_data);

  num_dofs=bcm_CSRMatrixNumRows(A);

  double time1, time2;
  *timetot=0.0;

  /* Create and Initialize Workspaces */
  d_vector= (bcm_Vector **) calloc(2, sizeof(bcm_Vector));
  for(i=0; i<2; i++) 
    {
      d_vector[i]=bcm_VectorCreate(num_dofs);
      ierr=bcm_VectorInitialize(d_vector[i]);
    }

  v_vector=bcm_VectorCreate(num_dofs);
  ierr=bcm_VectorInitialize(v_vector);
  w_vector=bcm_VectorCreate(num_dofs);
  ierr=bcm_VectorInitialize(w_vector);

  /* sparse-matrix vector product: --------------------------*/

  ierr=bcm_CSRMatrixMatvec(one,A,x_vector,zero,v_vector);

  /* compute initial residual: w <-- rhs - A *x = rhs - v;          */
  bcm_VectorCopy(rhs_vector,w_vector);
  ierr= bcm_VectorAxpy(minusone, v_vector, w_vector);

  delta0 = bcm_VectorNorm(w_vector);
  rhs_norm = bcm_VectorNorm(rhs_vector);

  if (delta0 <= eps * rhs_norm)
    {
      *num_iter = 0;
      return ierr;
    }

  idx=0;
  if(precon)
    {
      /* apply preconditioner to w */
      time1=time_getWallclockSeconds();
      ierr=bcm_PrecApply(bootamg_data, boot_amg, amg_cycle, w_vector, d_vector[idx]);
      time2=time_getWallclockSeconds()-time1;
      *timetot=*timetot+time2;
    }

  delta_old = bcm_VectorInnerProd(w_vector, d_vector[idx]);
  if (delta_old < zero)
    {
      printf("\n ERROR1: indefinite preconditioner in cg_iter_coarse: %e\n", delta_old);
      return -1;
    }

 loop:

  /* sparse-matrix vector product: --------------------------*/
      //time1=time_getWallclockSeconds();
      ierr=bcm_CSRMatrixMatvec(one,A,d_vector[idx],zero,v_vector);
      //time2=time_getWallclockSeconds()-time1;
      //printf("timematvec=%e\n",time2);

  tau = bcm_VectorInnerProd(d_vector[idx], v_vector);

  if (tau <= zero)
    {
      printf("\n ERROR2: indefinite matrix in cg_iter_coarse: %e\n", tau);
      return -1;
    }

  alpha = delta_old/tau;

  /* update solution  */
  ierr=bcm_VectorAxpy(alpha,d_vector[idx],x_vector);

  /* update residual  */
  ierr=bcm_VectorAxpy(-alpha,v_vector,w_vector);

  l2_norm = bcm_VectorNorm(w_vector);

  iter++;
  idx=iter%2;
  
  bcm_VectorSetConstantValues(d_vector[idx], zero);
  if(precon)
    {
      /* apply preconditioner to w */
      time1=time_getWallclockSeconds();
      ierr=bcm_PrecApply(bootamg_data, boot_amg, amg_cycle, w_vector, d_vector[idx]);
      time2=time_getWallclockSeconds()-time1;
      *timetot=*timetot+time2;
    }

  /* update direction  */
      tau1 = bcm_VectorInnerProd(d_vector[idx], v_vector);
      beta= tau1/tau;
      if(idx==1) ierr=bcm_VectorAxpy(-beta,d_vector[idx-1],d_vector[idx]);
      else ierr=bcm_VectorAxpy(-beta,d_vector[idx+1],d_vector[idx]);

      delta_old = bcm_VectorInnerProd(w_vector,d_vector[idx]);

  printf("bootpcg iteration: %d;  residual: %e, relative residual: %e\n", iter, l2_norm, l2_norm/delta0);
  if (l2_norm > rtol * delta0 && iter < max_iter) goto loop;
 
  *num_iter = iter;

  bcm_VectorDestroy(v_vector);
  bcm_VectorDestroy(w_vector);  
  for(i=0; i<2; i++) bcm_VectorDestroy(d_vector[i]);  
  free(d_vector);
  
  return ierr;
}
