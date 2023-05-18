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

#include "bcm_bootamg.h"


#ifdef HAVE_SUPERLU
#include "slu_ddefs.h"
#endif

/*--------------------------------------------------------------------------
 * bcm_GAMGCycle: this function applies a general cycle (V-W-H-K), depending
 * on the user choice to each hierarchy generated in the bootstrap AMG process
 *--------------------------------------------------------------------------*/

int bcm_GAMGCycle(int k, bcm_BootAMGBuildData *bootamg_data,
		  bcm_BootAMG *boot_amg, bcm_AMGApplyData *amg_cycle, bcm_Vector **Rhs,
		  bcm_Vector **Xtent, int l)
{

  /* Data Structure variables */

  bcm_AMGHierarchy **Harray;
  Harray=bcm_BootAMGHarray(boot_amg);

  bcm_CSRMatrix    **A_array;
  bcm_CSRMatrix    **P_array;
  bcm_CSRMatrix    **L_array;
  bcm_CSRMatrix    **U_array;
  bcm_Vector       **D_array;
  bcm_Vector       *Res;

  int        cycle_type, relax_type, coarse_solver;
  int        prerelax_number, postrelax_number;
  double     relax_weight;
  int        *num_grid_sweeps;

  double time1, time2,tnrm;
  int debug=0;
  /* Acquire data and allocate memory */

  A_array           = bcm_AMGHierarchyAArray(Harray[k]);
  P_array           = bcm_AMGHierarchyPArray(Harray[k]);
  L_array           = bcm_AMGHierarchyLArray(Harray[k]);
  U_array           = bcm_AMGHierarchyUArray(Harray[k]);
  D_array           = bcm_AMGHierarchyDArray(Harray[k]);
  int num_levels    = bcm_AMGHierarchyNumLevels(Harray[k]);
#ifdef HAVE_SUPERLU
  factors_t        *LUfactors=bcm_AMGHierarchyLUfactors(Harray[k]);
#endif
  

  bcm_AMGBuildData *amg_data;
  amg_data=bcm_BootAMGBuildDataCompData(bootamg_data);
  coarse_solver     = bcm_AMGBuildDataCoarseSolver(amg_data);

  relax_type        = bcm_AMGApplyDataRelaxType(amg_cycle);
  relax_weight      = bcm_AMGApplyDataRelaxWeight(amg_cycle); 
  prerelax_number   = bcm_AMGApplyDataPreRelax(amg_cycle);
  postrelax_number  = bcm_AMGApplyDataPostRelax(amg_cycle);
  cycle_type        = bcm_AMGApplyDataCycleType(amg_cycle);
  num_grid_sweeps   = bcm_AMGApplyDataGridSweeps(amg_cycle);

  /* Local variables  */
  int   i, ierr, sizeprol, sizeres;
  int istop;
  double rtol=0.25, eps = DBL_EPSILON;

  ierr=0;
  Res=bcm_VectorCreate(bcm_VectorSize(Rhs[l-1]));
  bcm_VectorInitialize(Res);

  int prerelaxtype, postrelaxtype;
  if(relax_type == 1) 
    {
      prerelaxtype=1;
      postrelaxtype=2;
    }
  else
    {
      prerelaxtype=relax_type;
      postrelaxtype=relax_type;
    }

  if (debug) fprintf(stderr,"GAMGCycle: start of level %d\n",l);
  if(l == num_levels) 
    {
#ifdef HAVE_SUPERLU
       if(coarse_solver == 9)
	{
	  /* use SuperLU */
	  /* Is it first time around ? */
	  SuperMatrix *L, *U, B;
	  int *perm_r; /* row permutations from partial pivoting */
	  int *perm_c; /* column permutation vector */
	  double *b;
	  int      i, nrhs=1, info;
	  trans_t   trans = NOTRANS;
	  SuperLUStat_t stat;	      
	  int  n    = bcm_CSRMatrixNumRows(A_array[l-1]);
	  int  nnz  = bcm_CSRMatrixNumNonzeros(A_array[l-1]);
	  if (LUfactors==NULL)  {
	    fprintf(stderr,
		    "Fatal error: LUfactors should have been initialized\n");
	    exit(-1);	      
	  }
	  bcm_VectorCopy(Rhs[l-1],Xtent[l-1]);
	  b = bcm_VectorData(Xtent[l-1]);
	  StatInit(&stat);
	  
	  /* Extract the LU factors in the factors handle */
	  L      = LUfactors->L;
	  U      = LUfactors->U;
	  perm_c = LUfactors->perm_c;
	  perm_r = LUfactors->perm_r;
	  dCreate_Dense_Matrix(&B, n, nrhs, b, n, SLU_DN, SLU_D, SLU_GE);
	  /* Solve the system A*X=B, overwriting B with X. */
	  dgstrs(trans, L, U, perm_c, perm_r, &B, &stat, &info);
	  
	  Destroy_SuperMatrix_Store(&B);
	  StatFree(&stat);
	    
	}
      else
#endif
	
	{
	  bcm_CSRMatrixRelax(A_array[l-1],L_array[l-1],U_array[l-1],D_array[l-1],
			     Rhs[l-1],coarse_solver,relax_weight,Xtent[l-1]);
	}
    }
  else 
    {
      /* presmoothing steps */
      for(i=1; i<= prerelax_number; i++)
	bcm_CSRMatrixRelax(A_array[l-1],L_array[l-1],U_array[l-1],D_array[l-1],
			   Rhs[l-1],prerelaxtype,relax_weight,Xtent[l-1]);
      if (debug) {
	tnrm = bcm_VectorNorm(Rhs[l-1]);
	fprintf(stderr,"RHS at level %d %lg\n",l,tnrm);
	tnrm = bcm_VectorNorm(Xtent[l-1]);
	fprintf(stderr,"After first smoother at level %d XTent %lg \n",l,tnrm);
      }
  
      /* compute residual */
      
      bcm_VectorCopy(Rhs[l-1],Res);
      bcm_CSRMatrixMatvec(-1.0,A_array[l-1],Xtent[l-1],1.0,Res);
      if (debug) {
	tnrm = bcm_VectorNorm(Res);
	fprintf(stderr,"Residual  at level %d %lg\n",l,tnrm);
      }
      
      /* restrict residual */
      bcm_CSRMatrixMatvecT(1.0,P_array[l-1],Res,0.0,Rhs[l]);
      if (debug) {
	tnrm = bcm_VectorNorm(Rhs[l]);
	fprintf(stderr,"Next RHS %d %lg\n",l+1,tnrm);
      }
      
      bcm_VectorSetConstantValues(Xtent[l],0.0);
      /* coarse correction */
      if(cycle_type==3)
	{
	  if(l==num_levels-1 )
	    {
	      bcm_GAMGCycle(k, bootamg_data, boot_amg, amg_cycle, Rhs, Xtent, l+1);
	    }
	  else
	    { 
	      time1=time_getWallclockSeconds();
	      istop= bcm_inneritkcycle(k, Xtent[l], Rhs[l], bootamg_data,
				       boot_amg,amg_cycle, rtol,l);    
	      time2=time_getWallclockSeconds()-time1;
	    }
	}
      else
	{
	  for(i=1; i<=num_grid_sweeps[l-1]; i++)
	    {
	      bcm_GAMGCycle(k, bootamg_data, boot_amg, amg_cycle, Rhs, Xtent, l+1);
	      if(l == num_levels-1) break;
	    }
	}
      /* prolongate error */
      bcm_CSRMatrixMatvec(1.0,P_array[l-1],Xtent[l],1.0,Xtent[l-1]);
      if (debug) {
	tnrm = bcm_VectorNorm(Xtent[l-1]);
	fprintf(stderr,"After recursion at level %d XTent %lg \n",l,tnrm);
      }
      /* postsmoothing steps */
      for(i=1; i<= postrelax_number; i++)
	bcm_CSRMatrixRelax(A_array[l-1],L_array[l-1],U_array[l-1],D_array[l-1],Rhs[l-1],
			   postrelaxtype,relax_weight,Xtent[l-1]);
      if (debug) {
	tnrm = bcm_VectorNorm(Xtent[l-1]);
	fprintf(stderr,"After second smoother at level %d XTent %lg \n",l,tnrm);
      }
    }
  
  bcm_VectorDestroy(Res);
  if (debug) fprintf(stderr,"GAMGCycle: end of level %d\n",l);
  return ierr;
}
