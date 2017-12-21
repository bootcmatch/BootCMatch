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

#include "bcm_bootamg.h"

/* This function applies 2 iterations of the Preconditioned FCG at each level
 * but the coarsest one
 * of a MG hierarchy; it is called by the bcm_GAMGCycle function. We considered
 * the algorithm described in Y. Notay "An Aggregation-based Algebraic Multigrid
 * Mehtod" (Algorithm 3.2), ETNA 2010, Vol. 37*/

int
bcm_inneritkcycle(int kh, bcm_Vector *x_vector, bcm_Vector *rhs_vector,
		  bcm_BootAMGBuildData *bootamg_data,
		  bcm_BootAMG * boot_amg, bcm_AMGApplyData *amg_cycle, double rtol, int l)
{

  int ierr=0, i, k;
  double rhs_norm ,delta0, delta_old, delta, eps = DBL_EPSILON; 
  bcm_Vector **d_vector, *v_vector, *w_vector, *v1_vector; 
  int num_dofs;
  /* flexible 1 CG */
  int iter=0, idx, kk; 
  double tau, tau1, alpha, tau2, tau3, tau4;
  int debug=0;
  double l2_norm,tnrm;
  double   zero=0.0, one=1.0, minusone=-1.0;

  /* matrix data */
  bcm_AMGHierarchy **Harray;
  Harray=bcm_BootAMGHarray(boot_amg);

  bcm_CSRMatrix **A;
  A = bcm_AMGHierarchyAArray(Harray[kh]);

  num_dofs=bcm_CSRMatrixNumRows(A[l]);

  if (debug) fprintf(stderr,"Start inneritkcyle level %d\n",l+1);
  /* Create and Initialize Workspaces */
  d_vector= (bcm_Vector **) calloc(2, sizeof(bcm_Vector));
  for(i=0; i<2; i++) 
    {
      d_vector[i]=bcm_VectorCreate(num_dofs);
      ierr=bcm_VectorInitialize(d_vector[i]);
    }

  v_vector=bcm_VectorCreate(num_dofs);
  ierr=bcm_VectorInitialize(v_vector);
  v1_vector=bcm_VectorCreate(num_dofs);
  ierr=bcm_VectorInitialize(v1_vector);
  w_vector=bcm_VectorCreate(num_dofs);
  ierr=bcm_VectorInitialize(w_vector);

  bcm_VectorCopy(rhs_vector,w_vector);
  delta0 = bcm_VectorNorm(w_vector);
  if (debug) fprintf(stderr,"level %d delta0  %g \n",l+1,delta0);
  idx=0;
  /* apply preconditioner to w */
  ierr=bcm_KCApply(kh,bootamg_data, boot_amg, amg_cycle, l+1, w_vector, d_vector[idx]);
  delta_old = bcm_VectorInnerProd(w_vector, d_vector[idx]);
  
  if (debug) {
    tnrm = bcm_VectorNorm(w_vector);
    fprintf(stderr,"level %d recursion output W nrm   %g \n",l+1,tnrm);
    tnrm = bcm_VectorNorm(d_vector[idx]);
    fprintf(stderr,"level %d recursion output nrm   %g \n",l+1,tnrm);
  }
  if (delta_old < zero)
    {
      printf("\n ERROR1: indefinite preconditioner in inner_iter: %e\n", delta_old);
      return -1;
    }

  /* sparse-matrix vector product: --------------------------*/
  ierr=bcm_CSRMatrixMatvec(one,A[l],d_vector[idx],zero,v_vector);

  tau = bcm_VectorInnerProd(d_vector[idx], v_vector);
  if (debug) fprintf(stderr,"level %d delta_old %g tau %g\n",l+1,delta_old,tau);

  if (tau <= zero)
    {
      printf("\n ERROR2: indefinite matrix in inner_iter: %e\n", tau);
      return -1;
    }

  alpha = delta_old/tau;
  /* update residual  */
  ierr=bcm_VectorAxpy(-alpha,v_vector,w_vector);
  
  l2_norm = bcm_VectorNorm(w_vector);
  if (debug) fprintf(stderr,"level %d alpha %g l2_n %g rtol*delta0 %g \n",
		     l+1,alpha,l2_norm,rtol*delta0);
  if (l2_norm <= rtol * delta0)
    { 
      /* update solution  */
      ierr=bcm_VectorAxpy(alpha,d_vector[idx],x_vector);
    }
  else
    { 
      iter++;
      idx=iter%2;
  
      /* apply preconditioner to w */
      ierr=bcm_KCApply(kh,bootamg_data, boot_amg,amg_cycle,l+1,w_vector, d_vector[idx]);

      /* sparse-matrix vector product: --------------------------*/
      ierr=bcm_CSRMatrixMatvec(one,A[l],d_vector[idx],zero,v1_vector);

      tau1 = bcm_VectorInnerProd(d_vector[idx], v_vector); /* gamma of Notay algorithm */
      tau2 = bcm_VectorInnerProd(d_vector[idx], v1_vector);/* beta of Notay algorithm */
      tau3 = bcm_VectorInnerProd(d_vector[idx],w_vector); /* alpha 2 of Notay algorithm */
      tau4 = tau2 - pow(tau1,2)/tau; /* rho2 of Notay algorihtm */
      if (debug) fprintf(stderr,"tau 1:4 %g %g %g %g \n",tau1,tau2,tau3,tau4);

      /* update solution  */
      alpha=alpha-(tau1*tau3)/(tau*tau4);
      ierr=bcm_VectorAxpy(alpha,d_vector[idx-1],x_vector);
      alpha=tau3/tau4;
      ierr=bcm_VectorAxpy(alpha,d_vector[idx],x_vector);
    }

  bcm_VectorDestroy(v_vector);
  bcm_VectorDestroy(v1_vector);
  bcm_VectorDestroy(w_vector);  
  for(i=0; i<2; i++) bcm_VectorDestroy(d_vector[i]);  
  free(d_vector);
  if (debug) fprintf(stderr,"End inneritkcyle level %d\n",l);  
  return ierr;
}
