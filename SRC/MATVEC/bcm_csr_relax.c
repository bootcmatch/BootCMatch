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
/******************************************************************************
 *
 * Relaxation (Smoother) schemes
 *
 *****************************************************************************/

#include "bcm_matvec.h"
#ifdef HAVE_SUPERLU
#include "slu_ddefs.h"
#endif
/*--------------------------------------------------------------------------
 * bcm_CSRMatrixRelax
 *--------------------------------------------------------------------------*/

int  bcm_CSRMatrixRelax( bcm_CSRMatrix *A, bcm_CSRMatrix *L, bcm_CSRMatrix *U, bcm_Vector *D,
                         bcm_Vector     *f,
                         int             relax_type,
                         double          relax_weight,
                         bcm_Vector    *u )
{
  int             n       = bcm_CSRMatrixNumRows(A);
  double zero=0.0, one=1.0;
  int             relax_error = 0;
   
  bcm_Vector *utemp, *w;
  int ierr;

  utemp=bcm_VectorCreate(n);
  bcm_VectorInitialize(utemp);
  w=bcm_VectorCreate(n);
  bcm_VectorInitialize(w);
 
  /*-----------------------------------------------------------------------
   *     Switch statement to direct control based on relax_type:
   *     relax_type = 0 -> Jacobi
   *     relax_type = 1 -> Gauss-Seidel forward
   *     relax_type = 2 -> Gauss-Seidel backward
   *     relax_type = 3 -> symm. Gauss-Seidel
   *     relax_type = 9 -> Direct Solve
   * -----------------------------------------------------------------------*/

   

  switch (relax_type)
    {
    case 0: /* Weighted Jacobi */
      {
      	ierr=bcm_CSRMatrixMatvec(one,A,u,zero,utemp);
     	bcm_VectorCopy(f,w);
     	ierr= bcm_VectorAxpy(-one, utemp, w);
	bcm_VectorDestroy(utemp);
        utemp=bcm_VectorDiagScal(w,D);
        ierr=bcm_VectorAxpy(relax_weight,utemp,u); 
      }
      break;
    case 1: /* Gauss-Seidel forward */
      {
        ierr=bcm_CSRMatrixMatvec(one,U,u,zero,utemp);
        bcm_VectorCopy(f,w);
        ierr= bcm_VectorAxpy(-one, utemp, w);
	bcm_VectorDestroy(utemp);
        utemp=bcm_CSRMatrixLSolve(L,D,w);
        bcm_VectorCopy(utemp,u);
      }
      break;
    case 2: /* Gauss-Seidel backward */
      {
	ierr=bcm_CSRMatrixMatvec(one,L,u,zero,utemp);
	bcm_VectorCopy(f,w);
	ierr= bcm_VectorAxpy(-one, utemp, w);
	bcm_VectorDestroy(utemp);
	utemp=bcm_CSRMatrixUSolve(U,D,w);
        bcm_VectorCopy(utemp,u);
      }
      break;
    case 3: /* symm. Gauss-Seidel */
      {
	ierr=bcm_CSRMatrixMatvec(one,U,u,zero,utemp);
	bcm_VectorCopy(f,w);
	ierr= bcm_VectorAxpy(-one, utemp, w);
	bcm_VectorDestroy(utemp);
        utemp=bcm_CSRMatrixLSolve(L,D,w);
        bcm_VectorCopy(utemp,u);

	ierr=bcm_CSRMatrixMatvec(one,L,u,zero,utemp);
	bcm_VectorCopy(f,w);
	ierr= bcm_VectorAxpy(-one, utemp, w);
	bcm_VectorDestroy(utemp);		
	utemp=bcm_CSRMatrixUSolve(U,D,w);
        bcm_VectorCopy(utemp,u);
      }
      break;

    case 9: /* Direct solve: use gaussian elimination */
      {

	int gselim(double *A_mat, double *b_vec, int n);

        double *A_mat;
        double *b_vec;

	int             i, ii;
	int             jj;
	int             column;

	double         *A_data  = bcm_CSRMatrixData(A);
	int            *A_i     = bcm_CSRMatrixI(A);
	int            *A_j     = bcm_CSRMatrixJ(A);

	int             n       = bcm_CSRMatrixNumRows(A);
	int             nnz     = bcm_CSRMatrixNumNonzeros(A);

	double         *u_data  = bcm_VectorData(u);
	double         *f_data  = bcm_VectorData(f);


#ifdef HAVE_SUPERLU

	bcm_VectorCopy(f,u);

	double         *values;
	int            *rowind,*colptr;
	SuperMatrix SA, AC,B;
	SuperMatrix *SL, *SU;
	int *perm_r; /* row permutations from partial pivoting */
	int *perm_c; /* column permutation vector */
	int *etree;  /* column elimination tree */
	SCformat *Lstore;
	NCformat *Ustore;
	int      panel_size, permc_spec, relax;
	trans_t  trans;
	mem_usage_t   mem_usage;
	superlu_options_t options;
	SuperLUStat_t stat;
	//factors_t *LUfactors;
	GlobalLU_t Glu;   /* Not needed on return. */
	int info;
	dCompRow_to_CompCol(n,  n, nnz,
			    A_data, A_j, A_i,
			    &values,&rowind,&colptr);
	trans = NOTRANS;
	
	
	/* Set the default input options. */
	set_default_options(&options);
	
	/* Initialize the statistics variables. */
	StatInit(&stat);
	
	dCreate_CompCol_Matrix(&SA, n, n, nnz, values, rowind, colptr,
			       SLU_NC, SLU_D, SLU_GE);
	SL = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
	SU = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
	if ( !(perm_r = intMalloc(n)) ) ABORT("Malloc fails for perm_r[].");
	if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
	if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
    /*
	 * Get column permutation vector perm_c[], according to permc_spec:
	 *   permc_spec = 0: natural ordering 
	 *   permc_spec = 1: minimum degree on structure of A'*A
	 *   permc_spec = 2: minimum degree on structure of A'+A
	 *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
	 */    	
	options.ColPerm=2;
	permc_spec = options.ColPerm;
	get_perm_c(permc_spec, &SA, perm_c);
	
	sp_preorder(&options, &SA, perm_c, etree, &AC);
	
	panel_size = sp_ienv(1);
	relax = sp_ienv(2);
#if defined(SLU_VERSION_5)
	dgstrf(&options, &AC, relax, panel_size, etree,
	       NULL, 0, perm_c, perm_r, SL, SU, &Glu, &stat, &info);
#elif defined(SLU_VERSION_4)
	dgstrf(&options, &AC, relax, panel_size, etree,
	       NULL, 0, perm_c, perm_r, SL, SU, &stat, &info);
#else
    choke_on_me;
#endif
    
    if ( info == 0 ) {
      Lstore = (SCformat *) SL->Store;
      Ustore = (NCformat *) SU->Store;
      dQuerySpace(SL, SU, &mem_usage);
#if 0
      printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
      printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
      printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz);
      printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
	     mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
#endif
    } else {
      printf("dgstrf() error returns INFO= %d\n", info);
      if ( info <= n ) { /* factorization completes */
	dQuerySpace(SL, SU, &mem_usage);
	printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
	       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
      }
    }
    dCreate_Dense_Matrix(&B, n, 1, u_data, n, SLU_DN, SLU_D, SLU_GE);
    /* Solve the system A*X=B, overwriting B with X. */
    dgstrs(trans, SL, SU, perm_c, perm_r, &B, &stat, &info);
    
    Destroy_SuperMatrix_Store(&B);
    SUPERLU_FREE(etree);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_SuperNode_Matrix(SL);
    Destroy_CompCol_Matrix(SU);
    Destroy_CompCol_Matrix(&SA);
    Destroy_CompCol_Permuted(&AC);
    SUPERLU_FREE (SL);
    SUPERLU_FREE (SU);


    StatFree(&stat);

	
#else
	A_mat = (double *) calloc(n*n, sizeof(double));
	b_vec = (double *) calloc(n, sizeof(double));    

	/*-----------------------------------------------------------------
	 *  Load CSR matrix into A_mat.
	 *-----------------------------------------------------------------*/

	for (i = 0; i < n; i++)
	  {
            for (jj = A_i[i]; jj < A_i[i+1]; jj++)
	      {
		column = A_j[jj];
		A_mat[i*n+column] = A_data[jj];
	      }
            b_vec[i] = f_data[i];
	  }

	relax_error = gselim(A_mat,b_vec,n);

	for (i = 0; i < n; i++)
	  {
            u_data[i] = b_vec[i];
	  }

	free(A_mat); 
	free(b_vec);
         
#endif
      }
      break;   
    }

  bcm_VectorDestroy(utemp);
  bcm_VectorDestroy(w);
  return(relax_error); 
}

/*-------------------------------------------------------------------------
 *
 *                      Gaussian Elimination
 *
 *------------------------------------------------------------------------ */

int gselim(A,x,n)
     double *A;
     double *x;
     int n;
{
  int    err_flag = 0;
  int    j,k,m;
  double factor, maxpivot;

  if (n==1)                           /* A is 1x1 */  
    {
      if (fabs(A[0]) >= DBL_EPSILON)
	{
	  x[0] = x[0]/A[0];
	  return(err_flag);
	}
      else
	{
	  err_flag = 1;
	  printf("Warning: null diagonal element in Gaussian Elimination\n");
	  return(err_flag);
	}
    }
  else                               /* A is nxn.  Forward elimination */ 
    {
      int *piv;
      piv = (int *) calloc(n, sizeof(int));
      for(k=0; k<n; k++) piv[k]=k;

      for (k = 0; k < n-1; k++)
	{

	  maxpivot=fabs(A[piv[k]*n+k]);
	  for (j = k+1; j < n; j++)
	    {
	      /* apply partial pivoting */
	      int tmp;
	      if( fabs(A[piv[j]*n+k]) > maxpivot) 
		{
		  tmp=piv[k];
		  piv[k]=piv[j];
		  piv[j]=tmp;
		  maxpivot=fabs(A[piv[j]*n+k]);
		}
	      if (fabs(A[piv[j]*n+k]) >= DBL_EPSILON)
		{
		  factor = A[piv[j]*n+k]/A[piv[k]*n+k];
		  for (m = k+1; m < n; m++)
                    {
		      A[piv[j]*n+m]  -= factor * A[piv[k]*n+m];
                    }
		  /* Elimination step for rhs */ 
		  x[piv[j]] -= factor * x[piv[k]];              
		}
	    }
	}
      /* Back Substitution  */
      for (k = n-1; k > 0; --k)
	{
	  if (fabs(A[piv[k]*n+k]) >= DBL_EPSILON) 
	    {
	      x[piv[k]] /= A[piv[k]*n+k];
	      for (j = 0; j < k; j++)
		{
                  if (fabs(A[piv[j]*n+k]) >= DBL_EPSILON) x[piv[j]] -= x[piv[k]] * A[piv[j]*n+k];
		}
	    }
	  else printf("Warning: singular matrix, null pivot in backward substitution, step %d\n", k);
	}
      if (fabs(A[piv[0]]) >= DBL_EPSILON) x[piv[0]] /= A[piv[0]];
      else printf("Warning: singular matrix, null pivot in backward substitution, step %d\n", 0);
      
      free(piv);
      return(err_flag);
    }
}
