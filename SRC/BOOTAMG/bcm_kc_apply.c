/*
                BootCMatch
     Bootstrap AMG based on Compatible weighted Matching, version 0.9
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

#include "bcm_bootamg.h"

/* ************************************************************************
*
*  function for applying a single AMG hierarchy (the kk component, starting 
*  from the level l), built within the bootstrap process, as preconditioner 
*  for a Krylov solver.
*
************************************************************************** */

int bcm_KCApply(int kk, bcm_BootAMGBuildData *bootamg_data,
		bcm_BootAMG *boot_amg, bcm_AMGApplyData *amg_cycle,
		int l, bcm_Vector *rhs, bcm_Vector *x)
{
   int ierr;
   int i, j, k;
   int num_lev;

   bcm_Vector **Xtent, **RHS;

   bcm_AMGHierarchy **Harray;
   Harray=bcm_BootAMGHarray(boot_amg);

   /* data related to single AMG */

   /* starting data */
   bcm_CSRMatrix **A_array;
   A_array = bcm_AMGHierarchyAArray(Harray[kk]);
   int nsize=bcm_CSRMatrixNumRows(A_array[l-1]);

/* Initialize data structure needed for cycling. 
   NB: Maybe these structures should be included in the amg_cycle structure! */
   num_lev=bcm_AMGHierarchyNumLevels(Harray[kk]);
   RHS = (bcm_Vector **) calloc(num_lev, sizeof(bcm_Vector));
   Xtent = (bcm_Vector **) calloc(num_lev, sizeof(bcm_Vector));
   
   for(i=l-1; i<num_lev; i++)
     {
       RHS[i]=bcm_VectorCreate(bcm_CSRMatrixNumRows(A_array[i]));
       Xtent[i]=bcm_VectorCreate(bcm_CSRMatrixNumRows(A_array[i]));
       bcm_VectorInitialize(RHS[i]);
       bcm_VectorInitialize(Xtent[i]);
     }
   bcm_VectorCopy(rhs,RHS[l-1]);
   
   bcm_GAMGCycle(kk,bootamg_data,boot_amg,amg_cycle, RHS, Xtent,l); 
   
   bcm_VectorCopy(Xtent[l-1],x);
   for(i=l-1; i<num_lev; i++)
     {
       bcm_VectorDestroy(RHS[i]);
       bcm_VectorDestroy(Xtent[i]);
     }
   free(RHS);
   free(Xtent);
   return ierr;
}
