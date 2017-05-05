/* 
                BootCMatch 
     Bootstrap AMG based on Compatible weighted Matching version 0.9
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
/******************************************************************************
 *
 * functions for AMG hierarchy data structures
 *
 *****************************************************************************/

#include "bcm_amg.h"
/*-----------------------------------------------------------------------
 * functions for bcm_AMGBuildData
 *-----------------------------------------------------------------------*/
void *
bcm_AMGInitialize()
{
   bcm_AMGBuildData *amg_data;

   /* setup params */
   int      maxlevels; /* max number of levels per hierarchy */
   int      maxcoarsesize; /* maximum size of the coarsest matrix */
   int      sweepnumber; /* number of pairwise aggregation steps */
   int      agg_interp_type; /* 1 for smoothed aggregation, 0 for pure aggregation */
   int      agg_match_type;  /* 1 for exact matching (HSL-mc64), 2 for auction based matching, 0 for Preis approximate matching */
   int      coarse_solver; /* solver to be used on the coarsest level */
   /*     relax_type/coarse_solver = 0 ->  1 sweep of Jacobi
    *     relax_type/coarse_solver = 1 ->  1 sweep of forward Gauss-Seidel
    *     relax_type/coarse_solver = 2 ->  1 sweep of backward Gauss-Seidel
    *     relax_type/coarse_solver = 3 ->  1 sweep of symm. Gauss-Seidel
    *     relax_type/coarse_solver = 9 -> Direct Solve */

   /* CR params */
   int      CRrelax_type; /* to choose relaxation scheme for Compatible Relaxation */
   int      CRit;  /* number of iterations for Compatible Relaxation */
   double   CRrelax_weight; /* weight for weighted Jacobi in CR */
   double   CRratio; /* optimal convergence ratio in Compatible Relaxation to stop coarsening*/

   /* Problem Data (to be read from user-defined function) */
   bcm_CSRMatrix *A; /* matrix */
   bcm_Vector *w;  /* arbitrary (smooth) vector */


   /*-----------------------------------------------------------------------
    * Setup default values for parameters
    *-----------------------------------------------------------------------*/

   /* setup params */
     maxlevels=100;
     maxcoarsesize=100; 
     sweepnumber=0; 
     agg_interp_type=0; 
     agg_match_type=0; 
     coarse_solver=9; 

   /* CR params */
     CRrelax_type=0; 
     CRit=20;  
     CRrelax_weight=1.0/3.0; 
     CRratio=0.3; 

   /*-----------------------------------------------------------------------
    * Create the bcm_AMGBuildData structure and return
    *-----------------------------------------------------------------------*/

   amg_data = (bcm_AMGBuildData *) calloc(1, sizeof(bcm_AMGBuildData));

   bcm_AMGBuildSetMaxLevels(amg_data, maxlevels);
   bcm_AMGBuildSetMaxCoarseSize(amg_data, maxcoarsesize);
   bcm_AMGBuildSetSweepNumber(amg_data, sweepnumber);
   bcm_AMGBuildSetAggInterpType(amg_data, agg_interp_type);
   bcm_AMGBuildSetAggMatchType(amg_data, agg_match_type);
   bcm_AMGBuildSetCoarseSolver(amg_data, coarse_solver);

   bcm_AMGBuildSetCRRelaxWeight(amg_data, CRrelax_weight);
   bcm_AMGBuildSetCRIterations(amg_data, CRit);
   bcm_AMGBuildSetCRRelaxType(amg_data, CRrelax_type);
   bcm_AMGBuildSetCRRatio(amg_data, CRratio );

   return (void *) amg_data;
}

/*--------------------------------------------------------------------------
 * Routines to set the setup phase parameters
 *--------------------------------------------------------------------------*/

int
bcm_AMGBuildSetMaxLevels( void *data,
                       int   max_levels )
{
   int ierr = 0;
   bcm_AMGBuildData  *amg_data = data;
 
   bcm_AMGBuildDataMaxLevels(amg_data) = max_levels;

   return (ierr);
}

int
bcm_AMGBuildSetAggInterpType( void     *data,
                        int       inter_type )
{
   int ierr = 0;
   bcm_AMGBuildData  *amg_data = data;

   bcm_AMGBuildDataAggInterpType(amg_data) = inter_type;

   return (ierr);
}

int
bcm_AMGBuildSetAggMatchType( void     *data,
                        int       match_type )
{
   int ierr = 0;
   bcm_AMGBuildData  *amg_data = data;

   bcm_AMGBuildDataAggMatchType(amg_data) = match_type;

   return (ierr);
}


int
bcm_AMGBuildSetMaxCoarseSize( void     *data,
                     int       maxcoarse_size )
{
   int ierr = 0;
   bcm_AMGBuildData  *amg_data = data;
 
   bcm_AMGBuildDataMaxCoarseSize(amg_data) = maxcoarse_size;

   return (ierr);
} 

int
bcm_AMGBuildSetSweepNumber( void     *data,
                           int      sweepnumber )
{
   int ierr = 0;
   bcm_AMGBuildData  *amg_data = data;

   bcm_AMGBuildDataSweepNumber(amg_data) = sweepnumber;

   return (ierr);
}
 
int
bcm_AMGBuildSetCRRelaxType( void     *data,
                           int      crrelax_type )
{
   int ierr = 0;
   bcm_AMGBuildData  *amg_data = data;

   bcm_AMGBuildDataCRRelaxType(amg_data) = crrelax_type;

   return (ierr);
}

int
bcm_AMGBuildSetCRRelaxWeight( void     *data,
                         double   crrelax_weight )
{
   int ierr = 0;
   bcm_AMGBuildData  *amg_data = data;

   bcm_AMGBuildDataCRRelaxWeight(amg_data) = crrelax_weight; 

   return (ierr);
}

int
bcm_AMGBuildSetCRIterations( void     *data,
                         int   criterations )
{
   int ierr = 0;
   bcm_AMGBuildData  *amg_data = data;

   bcm_AMGBuildDataCRIterations(amg_data) = criterations; 

   return (ierr);
}

int
bcm_AMGBuildSetCRRatio( void     *data,
                         double   crratio )
{
   int ierr = 0;
   bcm_AMGBuildData  *amg_data = data;

   bcm_AMGBuildDataCRratio(amg_data) = crratio; 

   return (ierr);
}

int
bcm_AMGBuildSetCoarseSolver( void     *data,
                     int      coarse_solver )
{
   int ierr = 0;
   bcm_AMGBuildData  *amg_data = data;
 
   bcm_AMGBuildDataCoarseSolver(amg_data) = coarse_solver;

   return (ierr);
} 

int 
bcm_AMGBuildDataDestroy(void *data)

{
   int ierr=0;
   bcm_AMGBuildData  *amg_data = data;

   if(amg_data)
   {
     free(amg_data);
   }
   return ierr;
}

/*-----------------------------------------------------------------------
 * functions for bcm_AMGHierarchy
 *-----------------------------------------------------------------------*/

bcm_AMGHierarchy *
bcm_AMGHierarchyCreate(int maxlevels)

{
  bcm_AMGHierarchy *amg_hierarchy;

  amg_hierarchy = (bcm_AMGHierarchy *) calloc(1, sizeof(bcm_AMGHierarchy));

  bcm_AMGHierarchyAArray(amg_hierarchy)=NULL;
  bcm_AMGHierarchyPArray(amg_hierarchy)=NULL;
  bcm_AMGHierarchyLArray(amg_hierarchy)=NULL;
  bcm_AMGHierarchyUArray(amg_hierarchy)=NULL;
  bcm_AMGHierarchyDArray(amg_hierarchy)=NULL;
#ifdef HAVE_SUPERLU
  bcm_AMGHierarchyLUfactors(amg_hierarchy)=NULL;
#endif

  bcm_AMGHierarchyNumLevels(amg_hierarchy)=maxlevels;
  bcm_AMGHierarchyOpCmplx(amg_hierarchy)=0;

  return amg_hierarchy;
}

int 
bcm_AMGHierarchyDestroy(bcm_AMGHierarchy *amg_hierarchy)

{
  int ierr=0, k, nl;
  bcm_CSRMatrix **tmpv;
  bcm_Vector **tmpv1;
#ifdef HAVE_SUPERLU
  factors_t *LUfactors;
#endif
  if(amg_hierarchy)
    {
      nl=bcm_AMGHierarchyNumLevels(amg_hierarchy);
      /* Do not destroy the fine matrix ! */
      tmpv =bcm_AMGHierarchyAArray(amg_hierarchy);
      for (k=1; k<nl; k++)
	bcm_CSRMatrixDestroy(tmpv[k]);
      free(tmpv);
      tmpv =bcm_AMGHierarchyLArray(amg_hierarchy);
      for (k=0; k<nl; k++)
	bcm_CSRMatrixDestroy(tmpv[k]);
      free(tmpv);
      tmpv =bcm_AMGHierarchyUArray(amg_hierarchy);
      for (k=0; k<nl; k++)
	bcm_CSRMatrixDestroy(tmpv[k]);
      free(tmpv);
      tmpv =bcm_AMGHierarchyPArray(amg_hierarchy);
      for (k=0; k<nl; k++)
	bcm_CSRMatrixDestroy(tmpv[k]);
      free(tmpv);
      tmpv1 =bcm_AMGHierarchyDArray(amg_hierarchy);
      for (k=0; k<nl; k++)
	bcm_VectorDestroy(tmpv1[k]);
      free(tmpv1);
#ifdef HAVE_SUPERLU
      LUfactors=bcm_AMGHierarchyLUfactors(amg_hierarchy);
      if (LUfactors != NULL) {
	SUPERLU_FREE (LUfactors->perm_r);
	SUPERLU_FREE (LUfactors->perm_c);
	Destroy_SuperNode_Matrix(LUfactors->L);
	Destroy_CompCol_Matrix(LUfactors->U);
	SUPERLU_FREE (LUfactors->L);
	SUPERLU_FREE (LUfactors->U);
	SUPERLU_FREE (LUfactors);
      }
#endif

       free(amg_hierarchy);
     }
   return ierr;
}

int 
bcm_AMGHierarchyInitialize(bcm_AMGHierarchy *amg_hierarchy)
{
   int num_levels=bcm_AMGHierarchyNumLevels(amg_hierarchy);

   int ierr=0;

   if(num_levels)
   {
   if(! bcm_AMGHierarchyAArray(amg_hierarchy))
      bcm_AMGHierarchyAArray(amg_hierarchy)= (bcm_CSRMatrix **) calloc(num_levels, sizeof(bcm_CSRMatrix));
   if(! bcm_AMGHierarchyPArray(amg_hierarchy))
      bcm_AMGHierarchyPArray(amg_hierarchy)= (bcm_CSRMatrix **) calloc(num_levels-1, sizeof(bcm_CSRMatrix));
   if(! bcm_AMGHierarchyLArray(amg_hierarchy))
      bcm_AMGHierarchyLArray(amg_hierarchy)= (bcm_CSRMatrix **) calloc(num_levels, sizeof(bcm_CSRMatrix));
   if(! bcm_AMGHierarchyUArray(amg_hierarchy))
      bcm_AMGHierarchyUArray(amg_hierarchy)= (bcm_CSRMatrix **) calloc(num_levels, sizeof(bcm_CSRMatrix));
   if(! bcm_AMGHierarchyDArray(amg_hierarchy))
      bcm_AMGHierarchyDArray(amg_hierarchy)= (bcm_Vector **) calloc(num_levels, sizeof(bcm_Vector));
   }
   return ierr;
}

/*-----------------------------------------------------------------------
 * functions for bcm_AMGApplyData
 *-----------------------------------------------------------------------*/

void *
bcm_AMGCycleInitialize()
{
   bcm_AMGApplyData *amg_cycle;

   /* cycle type params */
   int        cycle_type; /* 0 for V-cycle, 2 for W-cycle, 1 for G-cycle, 3 for k-cycle*/
   int        *num_grid_sweeps; /* number of sweeps on a fixed level in case of G-cycle */ 

   int        relax_type; /* type of pre/post relaxation/smoothing */
   int        prerelax_number; /* number of pre smoothing steps */
   int        postrelax_number; /* number of post smoothing steps */

   double     relax_weight; /* weight for Jacobi relaxation */

   /*-----------------------------------------------------------------------
    * Setup default values for parameters
    *-----------------------------------------------------------------------*/

   /* setup params */
     cycle_type=0;
     relax_type=0; 
     prerelax_number=1; 
     postrelax_number=1; 
     relax_weight=1.0; 
     num_grid_sweeps=NULL;

     
   /*-----------------------------------------------------------------------
    * Create the bcm_AMGApplyData structure and return
    *-----------------------------------------------------------------------*/

   amg_cycle = (bcm_AMGApplyData *) calloc(1, sizeof(bcm_AMGApplyData));

   bcm_AMGSetCycleType(amg_cycle, cycle_type);
   bcm_AMGSetNumGridSweeps( amg_cycle, num_grid_sweeps);
   bcm_AMGSetRelaxType(amg_cycle, relax_type);
   bcm_AMGSetPreRelaxSteps(amg_cycle, prerelax_number);
   bcm_AMGSetPostRelaxSteps(amg_cycle, postrelax_number);
   bcm_AMGSetRelaxWeight(amg_cycle, relax_weight);

   return (void *) amg_cycle;
}

/*--------------------------------------------------------------------------
 * Routines to set the setup phase parameters
 *--------------------------------------------------------------------------*/

int
bcm_AMGSetCycleType( void *data,
                       int   cycle_type )
{
   int ierr = 0;
   bcm_AMGApplyData  *amg_cycle = data;
 
   bcm_AMGApplyDataCycleType(amg_cycle) = cycle_type;

   return (ierr);
}

int
bcm_AMGSetNumGridSweeps( void     *data,
                        int      * num_grid_sweeps)
{
   int ierr = 0;
   bcm_AMGApplyData  *amg_cycle = data;
   
   bcm_AMGApplyDataGridSweeps(amg_cycle)=num_grid_sweeps;

   return (ierr);
}

int
bcm_AMGSetRelaxType( void     *data,
                     int      relax_type )
{
   int ierr = 0;
   bcm_AMGApplyData  *amg_cycle = data;
 
   bcm_AMGApplyDataRelaxType(amg_cycle) = relax_type;

   return (ierr);
} 


int
bcm_AMGSetPreRelaxSteps( void     *data,
                     int      prerelax_number )
{
   int ierr = 0;
   bcm_AMGApplyData  *amg_cycle = data;
 
   bcm_AMGApplyDataPreRelax(amg_cycle) = prerelax_number;

   return (ierr);
} 

int
bcm_AMGSetPostRelaxSteps( void     *data,
                     int      postrelax_number )
{
   int ierr = 0;
   bcm_AMGApplyData  *amg_cycle = data;
 
   bcm_AMGApplyDataPostRelax(amg_cycle) = postrelax_number;

   return (ierr);
} 

int
bcm_AMGSetRelaxWeight( void     *data,
                     double      relax_weight )
{
   int ierr = 0;
   bcm_AMGApplyData  *amg_cycle = data;
 
   bcm_AMGApplyDataRelaxWeight(amg_cycle) = relax_weight;

   return (ierr);
} 

int 
bcm_AMGApplyDataDestroy(void *data)

{
   int ierr=0;
   bcm_AMGApplyData  *amg_cycle = data;

   if(amg_cycle)
   {
     free(amg_cycle);
   }
   return ierr;
}
