/*
                BootCMatch
     Bootstrap AMG based on Compatible Matching version 1.0
    (C) Copyright 2017
                       Pasqua D'Ambra    ICAR-CNR
                       Salvatore Filippone Cranfield University
                       Panayot S. Vassilevski CACR-LLNL

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
*  function for applying a composite AMG, built by the bootstrap process,
*  as preconditioner of a Krylov solver
*
************************************************************************** */

int bcm_PrecApply(bcm_BootAMGBuildData *bootamg_data, bcm_BootAMG *boot_amg, bcm_AMGApplyData *amg_cycle, bcm_Vector *rhs, bcm_Vector *x)
{
  int ierr;
  double alpha;
  int n_hrc, solver_type;
  int i, k;
  int num_lev;
  
  bcm_Vector **Xtent, **RHS;
  
  bcm_AMGHierarchy **Harray;
  Harray=bcm_BootAMGHarray(boot_amg);
  n_hrc=bcm_BootAMGNHrc(boot_amg);
  
  /* data related to single AMG */
  bcm_CSRMatrix **A_array;
  A_array = bcm_AMGHierarchyAArray(Harray[0]);
  int nsize=bcm_CSRMatrixNumRows(A_array[0]);

  solver_type=bcm_BootAMGBuildDataCompType(bootamg_data);

  /* apply composite solver depending on the chosen type */
  switch(solver_type)
    {
    case 0: /* multiplicative */
      {
	for(k=0; k<n_hrc; k++)
	  {
	    /* Initialize cyclying structures */
	    num_lev=bcm_AMGHierarchyNumLevels(Harray[k]);
	    A_array=bcm_AMGHierarchyAArray(Harray[k]);
	    RHS = (bcm_Vector **) calloc(num_lev, sizeof(bcm_Vector));
	    Xtent = (bcm_Vector **) calloc(num_lev, sizeof(bcm_Vector));
	    for(i=0; i<num_lev; i++)
	      {
		RHS[i]=bcm_VectorCreate(bcm_CSRMatrixNumRows(A_array[i]));
		Xtent[i]=bcm_VectorCreate(bcm_CSRMatrixNumRows(A_array[i]));
		bcm_VectorInitialize(RHS[i]);
		bcm_VectorInitialize(Xtent[i]);
	      }
	    bcm_VectorCopy(rhs, RHS[0]);
	    bcm_VectorCopy(x, Xtent[0]);

	    bcm_GAMGCycle(k,bootamg_data,boot_amg,amg_cycle,RHS,Xtent,1);

	    bcm_VectorCopy(Xtent[0],x);

	    for(i=0; i<num_lev; i++)
	      {
		bcm_VectorDestroy(RHS[i]);
		bcm_VectorDestroy(Xtent[i]);
	      }
	    free(RHS);
	    free(Xtent);
	  }
      }
      break;
    case 1: /* symmetrized multiplicative */
      {
	for(k=0; k<n_hrc; k++)
	  {
	    /* Initialize cyclying structures */
	    num_lev=bcm_AMGHierarchyNumLevels(Harray[k]);
	    A_array=bcm_AMGHierarchyAArray(Harray[k]);
	    RHS = (bcm_Vector **) calloc(num_lev, sizeof(bcm_Vector));
	    Xtent = (bcm_Vector **) calloc(num_lev, sizeof(bcm_Vector));
	    for(i=0; i<num_lev; i++)
	      {
		RHS[i]=bcm_VectorCreate(bcm_CSRMatrixNumRows(A_array[i]));
		Xtent[i]=bcm_VectorCreate(bcm_CSRMatrixNumRows(A_array[i]));
		bcm_VectorInitialize(RHS[i]);
		bcm_VectorInitialize(Xtent[i]);
	      }
	    bcm_VectorCopy(rhs, RHS[0]);
	    bcm_VectorCopy(x, Xtent[0]);

	    bcm_GAMGCycle(k,bootamg_data,boot_amg,amg_cycle,RHS,Xtent,1);

	    bcm_VectorCopy(Xtent[0],x);

	    for(i=0; i<num_lev; i++)
	      {
		bcm_VectorDestroy(RHS[i]);
		bcm_VectorDestroy(Xtent[i]);
	      }
	    free(RHS);
	    free(Xtent);
	  }
	for(k=n_hrc-1; k>=0; k--)
	  {
	    /* Initialize cyclying structures */
	    num_lev=bcm_AMGHierarchyNumLevels(Harray[k]);
	    A_array=bcm_AMGHierarchyAArray(Harray[k]);
	    RHS = (bcm_Vector **) calloc(num_lev, sizeof(bcm_Vector));
	    Xtent = (bcm_Vector **) calloc(num_lev, sizeof(bcm_Vector));
	    for(i=0; i<num_lev; i++)
	      {
		RHS[i]=bcm_VectorCreate(bcm_CSRMatrixNumRows(A_array[i]));
		Xtent[i]=bcm_VectorCreate(bcm_CSRMatrixNumRows(A_array[i]));
		bcm_VectorInitialize(RHS[i]);
		bcm_VectorInitialize(Xtent[i]);
	      }
	    bcm_VectorCopy(rhs, RHS[0]);
	    bcm_VectorCopy(x, Xtent[0]);

	    bcm_GAMGCycle(k,bootamg_data,boot_amg,amg_cycle,RHS,Xtent,1);

	    bcm_VectorCopy(Xtent[0],x);
	    for(i=0; i<num_lev; i++)
	      {
		bcm_VectorDestroy(RHS[i]);
		bcm_VectorDestroy(Xtent[i]);
	      }
	    free(RHS);
	    free(Xtent);
	  }
      }
      break;
    case 2:  /* additive */
      {
	bcm_Vector *xadd;
	xadd=bcm_VectorCreate(nsize);
	bcm_VectorInitialize(xadd);
	bcm_VectorSetConstantValues(xadd,0.0);

	for(k=0; k<n_hrc; k++)
	  {
	    /* Initialize cyclying structures */
	    num_lev=bcm_AMGHierarchyNumLevels(Harray[k]);
	    A_array=bcm_AMGHierarchyAArray(Harray[k]);
	    RHS = (bcm_Vector **) calloc(num_lev, sizeof(bcm_Vector));
	    Xtent = (bcm_Vector **) calloc(num_lev, sizeof(bcm_Vector));
	    for(i=0; i<num_lev; i++)
	      {
		RHS[i]=bcm_VectorCreate(bcm_CSRMatrixNumRows(A_array[i]));
		Xtent[i]=bcm_VectorCreate(bcm_CSRMatrixNumRows(A_array[i]));
		bcm_VectorInitialize(RHS[i]);
		bcm_VectorInitialize(Xtent[i]);
	      }
	    bcm_VectorCopy(rhs, RHS[0]);
	    bcm_VectorCopy(x, Xtent[0]);

	    bcm_GAMGCycle(k,bootamg_data,boot_amg,amg_cycle,RHS,Xtent,1);

	    bcm_VectorAxpy(1.0,Xtent[0],xadd);
	    for(i=0; i<num_lev; i++)
	      {
		bcm_VectorDestroy(RHS[i]);
		bcm_VectorDestroy(Xtent[i]);
	      }
	    free(RHS);
	    free(Xtent);
	  }
	alpha=1.0/n_hrc;
	bcm_VectorScale(alpha,xadd);
	bcm_VectorCopy(xadd,x); 
	bcm_VectorDestroy(xadd); 
      }
      break; 
    }
     
  return ierr;
}
