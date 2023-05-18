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
/******************************************************************************
 *
 * functions for sparse matrix operations
 *
 *****************************************************************************/

#include "bcm_matvec.h"

/*--------------------------------------------------------------------------
 * bcm_CSRMatrixAdd:
 * adds two CSR Matrices A and B and returns a CSR Matrix C;
 * Note: The routine does not check for 0-elements which might be generated
 *       through cancellation of elements in A and B or already contained
 	 in A and B. To remove those, use bcm_CSRMatrixDeleteZeros 
 *--------------------------------------------------------------------------*/

bcm_CSRMatrix *
bcm_CSRMatrixAdd( bcm_CSRMatrix *A,
              bcm_CSRMatrix *B)
{
   double           *A_data   = bcm_CSRMatrixData(A);
   int              *A_i      = bcm_CSRMatrixI(A);
   int              *A_j      = bcm_CSRMatrixJ(A);
   int              nrows_A  =  bcm_CSRMatrixNumRows(A);
   int              ncols_A  =  bcm_CSRMatrixNumCols(A);
   double           *B_data   = bcm_CSRMatrixData(B);
   int              *B_i      = bcm_CSRMatrixI(B);
   int              *B_j      = bcm_CSRMatrixJ(B);
   int              nrows_B  = bcm_CSRMatrixNumRows(B);
   int              ncols_B  = bcm_CSRMatrixNumCols(B);
   bcm_CSRMatrix *C;
   double        *C_data;
   int	         *C_i;
   int           *C_j;

   int         ia, ib, ic, jcol, num_nonzeros;
   int	       pos;
   int         *marker;

   if (nrows_A != nrows_B || ncols_A != ncols_B)
     {
       printf("Warning! incompatible matrix dimensions!\n");
       return NULL;
     }


   marker = (int *) calloc(ncols_A, sizeof(int));
   C_i = (int *) calloc(nrows_A+1, sizeof(int));

   for (ia = 0; ia < ncols_A; ia++)
     marker[ia] = -1;

   num_nonzeros = 0;
   C_i[0] = 0;
   for (ic = 0; ic < nrows_A; ic++)
     {
       for (ia = A_i[ic]; ia < A_i[ic+1]; ia++)
	 {
	   jcol = A_j[ia];
	   marker[jcol] = ic;
	   num_nonzeros++;
	 }
       for (ib = B_i[ic]; ib < B_i[ic+1]; ib++)
	 {
	   jcol = B_j[ib];
	   if (marker[jcol] != ic)
	     {
	       marker[jcol] = ic;
	       num_nonzeros++;
	     }
	 }
       C_i[ic+1] = num_nonzeros;
     }

   C = bcm_CSRMatrixCreate(nrows_A, ncols_A, num_nonzeros);
   bcm_CSRMatrixI(C) = C_i;
   bcm_CSRMatrixInitialize(C);
   C_j = bcm_CSRMatrixJ(C);
   C_data = bcm_CSRMatrixData(C);

   for (ia = 0; ia < ncols_A; ia++)
     marker[ia] = -1;

   pos = 0;
   for (ic = 0; ic < nrows_A; ic++)
     {
       for (ia = A_i[ic]; ia < A_i[ic+1]; ia++)
	 {
	   jcol = A_j[ia];
	   C_j[pos] = jcol;
	   C_data[pos] = A_data[ia];
	   marker[jcol] = pos;
	   pos++;
	 }
       for (ib = B_i[ic]; ib < B_i[ic+1]; ib++)
	 {
	   jcol = B_j[ib];
	   if (marker[jcol] < C_i[ic])
	     {
	       C_j[pos] = jcol;
	       C_data[pos] = B_data[ib];
	       marker[jcol] = pos;
	       pos++;
	     }
	   else
	     {
	       C_data[marker[jcol]] += B_data[ib];
	     }
	 }
     }

   free(marker);
   return C;
}	

/*--------------------------------------------------------------------------
 * bcm_CSRMatrixMultiply
 * multiplies two CSR Matrices A and B and returns a CSR Matrix C;
 * Note: The routine does not check for 0-elements which might be generated
 *       through cancellation of elements in A and B or already contained
 	 in A and B. To remove those, use bcm_CSRMatrixDeleteZeros 
 *--------------------------------------------------------------------------*/

bcm_CSRMatrix *
bcm_CSRMatrixMultiply( bcm_CSRMatrix *A,
              bcm_CSRMatrix *B)
{
   double     *A_data   = bcm_CSRMatrixData(A);
   int        *A_i      = bcm_CSRMatrixI(A);
   int        *A_j      = bcm_CSRMatrixJ(A);
   int         nrows_A  = bcm_CSRMatrixNumRows(A);
   int         ncols_A  = bcm_CSRMatrixNumCols(A);
   double     *B_data   = bcm_CSRMatrixData(B);
   int        *B_i      = bcm_CSRMatrixI(B);
   int        *B_j      = bcm_CSRMatrixJ(B);
   int         nrows_B  = bcm_CSRMatrixNumRows(B);
   int         ncols_B  = bcm_CSRMatrixNumCols(B);
   bcm_CSRMatrix *C;
   double     *C_data;
   int	      *C_i;
   int        *C_j;

   int         ia, ib, ic, ja, jb, num_nonzeros=0;
   int	       row_start, counter;
   double      a_entry, b_entry;
   int         *B_marker;

   if (ncols_A != nrows_B)
     {
       printf("Warning! incompatible matrix dimensions!\n");
       return NULL;
     }


   B_marker = (int *) calloc(ncols_B, sizeof(int));
   C_i = (int *) calloc(nrows_A+1, sizeof(int));

   for (ib = 0; ib < ncols_B; ib++)
     B_marker[ib] = -1;

   for (ic = 0; ic < nrows_A; ic++)
     {
       for (ia = A_i[ic]; ia < A_i[ic+1]; ia++)
	 {
	   ja = A_j[ia];
	   for (ib = B_i[ja]; ib < B_i[ja+1]; ib++)
	     {
	       jb = B_j[ib];
	       if (B_marker[jb] != ic)
		 {
		   B_marker[jb] = ic;
		   num_nonzeros++;
		 }
	     }
	 }
       C_i[ic+1] = num_nonzeros;
     }

   C = bcm_CSRMatrixCreate(nrows_A, ncols_B, num_nonzeros);
   bcm_CSRMatrixI(C) = C_i;
   bcm_CSRMatrixInitialize(C);
   C_j = bcm_CSRMatrixJ(C);
   C_data = bcm_CSRMatrixData(C);

   for (ib = 0; ib < ncols_B; ib++)
     B_marker[ib] = -1;

   counter = 0;
   for (ic = 0; ic < nrows_A; ic++)
     {
       row_start = C_i[ic];
       for (ia = A_i[ic]; ia < A_i[ic+1]; ia++)
	 {
	   ja = A_j[ia];
	   a_entry = A_data[ia];
	   for (ib = B_i[ja]; ib < B_i[ja+1]; ib++)
	     {
	       jb = B_j[ib];
	       b_entry = B_data[ib];
	       if (B_marker[jb] < row_start)
		 {
		   B_marker[jb] = counter;
		   C_j[B_marker[jb]] = jb;
		   C_data[B_marker[jb]] = a_entry*b_entry;
		   counter++;
		 }
	       else
		 C_data[B_marker[jb]] += a_entry*b_entry;
				 
	     }
	 }
     }
   free(B_marker);
   bcm_CSRMatrixSort(C);
   return C;
}	

bcm_CSRMatrix *
bcm_CSRMatrixDeleteZeros( bcm_CSRMatrix *A, double tol)
{
  double     *A_data   = bcm_CSRMatrixData(A);
  int        *A_i      = bcm_CSRMatrixI(A);
  int        *A_j      = bcm_CSRMatrixJ(A);
  int         nrows_A  = bcm_CSRMatrixNumRows(A);
  int         ncols_A  = bcm_CSRMatrixNumCols(A);
  int         num_nonzeros  = bcm_CSRMatrixNumNonzeros(A);

  bcm_CSRMatrix *B;
  double     *B_data; 
  int        *B_i;
  int        *B_j;

  int zeros;
  int i, j;
  int pos_A, pos_B;

  zeros = 0;
  for (i=0; i < num_nonzeros; i++)
    if (fabs(A_data[i]) <= tol)
      zeros++;

  if (zeros)
    {
      B = bcm_CSRMatrixCreate(nrows_A,ncols_A,num_nonzeros-zeros);
      bcm_CSRMatrixInitialize(B);
      B_i = bcm_CSRMatrixI(B);
      B_j = bcm_CSRMatrixJ(B);
      B_data = bcm_CSRMatrixData(B);
      B_i[0] = 0;
      pos_A = 0;
      pos_B = 0;
      for (i=0; i < nrows_A; i++)
	{
	  for (j = A_i[i]; j < A_i[i+1]; j++)
	    {
	      if (fabs(A_data[j]) <= tol)
		{
		  pos_A++;
		}
	      else
		{
		  B_data[pos_B] = A_data[pos_A];
		  B_j[pos_B] = A_j[pos_A];
		  pos_B++;
		  pos_A++;
		}
	    }
	  B_i[i+1] = pos_B;
	}
      return B;
    }
  else
    return NULL;
}	


/******************************************************************************
 *
 * Compute transpose of a bcm_CSRMatrix
 *
 *****************************************************************************/

/*--------------------------------------------------------------------------
 * bcm_CSRMatrixTranspose
 *--------------------------------------------------------------------------*/

int bcm_CSRMatrixTranspose(bcm_CSRMatrix   *A, bcm_CSRMatrix   **AT,
			   int data)

{
  double       *A_data = bcm_CSRMatrixData(A);
  int          *A_i = bcm_CSRMatrixI(A);
  int          *A_j = bcm_CSRMatrixJ(A);
  int           num_rowsA = bcm_CSRMatrixNumRows(A);
  int           num_colsA = bcm_CSRMatrixNumCols(A);
  int           num_nonzerosA = bcm_CSRMatrixNumNonzeros(A);

  double       *AT_data;
  int          *AT_i;
  int          *AT_j;
  int           num_rowsAT;
  int           num_colsAT;
  int           num_nonzerosAT;

  int           max_col;
  int           i, j;

  /*-------------------------------------------------------------- 
   * First, ascertain that num_cols and num_nonzeros has been set. 
   * If not, set them.
   *--------------------------------------------------------------*/

  if (! num_nonzerosA)
    {
      num_nonzerosA = A_i[num_rowsA];
    }

  if (num_rowsA && ! num_colsA)
    {
      max_col = -1;
      for (i = 0; i < num_rowsA; ++i)
	{
          for (j = A_i[i]; j < A_i[i+1]; j++)
	    {
              if (A_j[j] > max_col)
		max_col = A_j[j];
	    }
	}
      num_colsA = max_col+1;
    }

  num_rowsAT = num_colsA;
  num_colsAT = num_rowsA;
  num_nonzerosAT = num_nonzerosA;

  *AT = bcm_CSRMatrixCreate(num_rowsAT, num_colsAT, num_nonzerosAT);

  AT_i = (int*) calloc(num_rowsAT+1, sizeof(int));
  AT_j = (int*) calloc(num_nonzerosAT, sizeof(int));
  bcm_CSRMatrixI(*AT) = AT_i;
  bcm_CSRMatrixJ(*AT) = AT_j;
  if (data) 
    {
      AT_data = (double *) calloc(num_nonzerosAT, sizeof(double));
      bcm_CSRMatrixData(*AT) = AT_data;
    }

  /*-----------------------------------------------------------------
   * Count the number of entries in each column of A (row of AT)
   * and fill the AT_i array.
   *-----------------------------------------------------------------*/

  for (i = 0; i < num_nonzerosA; i++)
    {
      ++AT_i[A_j[i]+1];
    }

  for (i = 2; i <= num_rowsAT; i++)
    {
      AT_i[i] += AT_i[i-1];
    }

  /*----------------------------------------------------------------
   * Load the data and column numbers of AT
   *----------------------------------------------------------------*/

  for (i = 0; i < num_rowsA; i++)
    {
      for (j = A_i[i]; j < A_i[i+1]; j++)
	{
	  assert(AT_i[A_j[j]] >= 0); 
	  assert(AT_i[A_j[j]] < num_nonzerosAT); 
	  AT_j[AT_i[A_j[j]]] = i;
	  if (data) AT_data[AT_i[A_j[j]]] = A_data[j];
	  AT_i[A_j[j]]++;
	}
    }

  /*------------------------------------------------------------
   * AT_i[j] now points to the *end* of the jth row of entries
   * instead of the beginning.  Restore AT_i to front of row.
   *------------------------------------------------------------*/

  for (i = num_rowsAT; i > 0; i--)
    {
      AT_i[i] = AT_i[i-1];
    }

  AT_i[0] = 0;

  return(0);
}

/*--------------------------------------------------------------------------
 * bcm_CSRMatrixSort:
 * Make sure the entries in each row are sorted in column-index order.
 *--------------------------------------------------------------------------*/
int bcm_CSRMatrixSort ( bcm_CSRMatrix *A )
{
  int i,j,nz;
  int *pja;
  double *pv;
  double       *A_data = bcm_CSRMatrixData(A);
  int          *A_i = bcm_CSRMatrixI(A);
  int          *A_j = bcm_CSRMatrixJ(A);
  int           num_rowsA = bcm_CSRMatrixNumRows(A);
  //return(0);
  for (i=0; i<num_rowsA; i++) {
    j   = A_i[i];
    nz  = A_i[i+1] - A_i[i];
    pja = &(A_j[j]);
    pv  = &(A_data[j]);
    bcm_SortMatrixRow(nz,pja,pv);
  }
  
  return(0);
}

/*--------------------------------------------------------------------------
 * bcm_SortMatrixRow:
 * Make sure the entries in a given row are sorted in column-index order.
 *--------------------------------------------------------------------------*/
int bcm_SortMatrixRow(int n, int *ja, double *dv)
{
  int i,j, jx;
  double dx;

  for (j=n-2; j>=0; j--) {
    if (ja[j+1] < ja[j]) {
      jx = ja[j];
      dx = dv[j];
      for (i=j+1; ; ) {
	ja[i-1]=ja[i];
	dv[i-1]=dv[i];
	i++;
	if (i>=n) break;
	if (ja[i] >=jx ) break;
      }
      ja[i-1] = jx;
      dv[i-1] = dx;
    }
  }  

  return(0); 
}


/*--------------------------------------------------------------------------
 * bcm_CSRMatrixReorder:
 * Reorders the column and data arrays of a square CSR matrix, such that the
 * first entry in each row is the diagonal one.
 *--------------------------------------------------------------------------*/

int bcm_CSRMatrixReorder(bcm_CSRMatrix *A)
{
  int i, j, tempi, row_size;
  double tempd;

  double *A_data = bcm_CSRMatrixData(A);
  int    *A_i = bcm_CSRMatrixI(A);
  int    *A_j = bcm_CSRMatrixJ(A);
  int     num_rowsA = bcm_CSRMatrixNumRows(A);
  int     num_colsA = bcm_CSRMatrixNumCols(A);

  /* the matrix should be square */
  if (num_rowsA != num_colsA)
    return -1;

  for (i = 0; i < num_rowsA; i++)
    {
      row_size = A_i[i+1]-A_i[i];

      for (j = 0; j < row_size; j++)
	{
	  if (A_j[j] == i)
	    {
	      if (j != 0)
		{
		  tempi = A_j[0];
		  A_j[0] = A_j[j];
		  A_j[j] = tempi;

		  tempd = A_data[0];
		  A_data[0] = A_data[j];
		  A_data[j] = tempd;
		}
	      break;
	    }

	  /* diagonal element is missing */
	  if (j == row_size-1)
            return -2;
	}

      A_j    += row_size;
      A_data += row_size;
    }

  return 0;
}

/*--------------------------------------------------------------------------
 * bcm_CSRMatrixSumElts:
 * Returns the sum of all matrix elements.
 *--------------------------------------------------------------------------*/

double bcm_CSRMatrixSumElts( bcm_CSRMatrix *A )
{
  double sum = 0;
  double * data = bcm_CSRMatrixData( A );
  int num_nonzeros = bcm_CSRMatrixNumNonzeros(A);
  int i;

  for ( i=0; i<num_nonzeros; ++i ) sum += data[i];

  return sum;
}

/*--------------------------------------------------------------------------
 * bcm_CSRMatrixInfNorm:
 * Returns the max norm of the matrix
 *--------------------------------------------------------------------------*/

double bcm_CSRMatrixInfNorm( bcm_CSRMatrix *A )
{
  double *sum, norm = 0.0;
  double * data = bcm_CSRMatrixData( A );
  int * A_i = bcm_CSRMatrixI( A );
  int * A_j = bcm_CSRMatrixJ( A );
  int num_rows = bcm_CSRMatrixNumRows(A);
  int i, j;

  sum=(double *) calloc(num_rows, sizeof(double));
  for ( i=0; i<num_rows; ++i ) 
    {
      for(j=A_i[i]; j< A_i[i+1]; j++) sum[i]=sum[i]+fabs(data[j]);
    }
  for ( i=0; i<num_rows; ++i ) 
    {
      if(sum[i] > norm) norm=sum[i];
    }
  free(sum);
  return norm;
}

/*--------------------------------------------------------------------------
 * bcm_VectorANorm:
 * Returns the A norm of the vector, where A is required to be a s.p.d. matrix.
 *--------------------------------------------------------------------------*/

double bcm_VectorANorm( bcm_CSRMatrix *A, bcm_Vector *vector )
{
  double norm = 0;
  bcm_Vector *temp;
   
  int n=bcm_VectorSize(vector);
  temp=bcm_VectorCreate(n);
  bcm_VectorInitialize(temp);

  bcm_CSRMatrixMatvec(1.0,A,vector,0.0,temp);
  norm=bcm_VectorInnerProd(temp,vector);
  norm=sqrt(norm);
  bcm_VectorDestroy(temp);
  return norm;
}

/*--------------------------------------------------------------------------
 * bcm_CSRMatrixLSolve(L,D,f):
 * Returns the solution of the lower triangular system (D+L)u=f 
 *--------------------------------------------------------------------------*/

bcm_Vector * bcm_CSRMatrixLSolve(bcm_CSRMatrix *L, bcm_Vector *D, bcm_Vector *f)
{
  int n = bcm_CSRMatrixNumRows(L);
  double *L_data=bcm_CSRMatrixData(L);
  double *D_data=bcm_VectorData(D);
  int *L_i=bcm_CSRMatrixI(L);
  int *L_j=bcm_CSRMatrixJ(L);

  double *f_data=bcm_VectorData(f);

  bcm_Vector *w;
  w=bcm_VectorCreate(n);
  bcm_VectorInitialize(w);

  double *w_data;
  w_data= bcm_VectorData(w);

  int i,j;
  int idx;

  w_data[0]=f_data[0]/D_data[0];

  for (i=1; i<n; ++i) 
    {
      for(j=L_i[i]; j< L_i[i+1]; j++) 
	{
	  idx=L_j[j];
	  w_data[i]=w_data[i]+L_data[j]*w_data[idx];
	}
      w_data[i]=(f_data[i]-w_data[i])/D_data[i];
    }
  return w;
}

/*--------------------------------------------------------------------------
 * bcm_CSRMatrixUSolve(U,D,f):
 * Returns the solution of the upper triangular system (U+D)u=f 
 *--------------------------------------------------------------------------*/

bcm_Vector * bcm_CSRMatrixUSolve(bcm_CSRMatrix *U, bcm_Vector *D, bcm_Vector *f)
{

  int n = bcm_CSRMatrixNumRows(U);
  double *U_data=bcm_CSRMatrixData(U);
  double *D_data=bcm_VectorData(D);
  int *U_i=bcm_CSRMatrixI(U);
  int *U_j=bcm_CSRMatrixJ(U);

  double *f_data=bcm_VectorData(f);

  bcm_Vector *w;
  w=bcm_VectorCreate(n);
  bcm_VectorInitialize(w);

  double *w_data;
  w_data= bcm_VectorData(w);

  int i,j,idx;

  w_data[n-1]=f_data[n-1]/D_data[n-1];

  for (i=n-2; i>=0; --i) 
    {
      for(j=U_i[i]; j<U_i[i+1]; j++)
	{
	  idx=U_j[j];
	  w_data[i]=w_data[i]+U_data[j]*w_data[idx];
	}
      w_data[i]=(f_data[i]-w_data[i])/D_data[i];
    }
  return w;
}
