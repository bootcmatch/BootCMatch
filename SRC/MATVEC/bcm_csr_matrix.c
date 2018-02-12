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
 * functions for bcm_CSRMatrix structure
 *
 *****************************************************************************/

#include "bcm_matvec.h"

/*--------------------------------------------------------------------------
 * bcm_CSRMatrixCreate
 *--------------------------------------------------------------------------*/

bcm_CSRMatrix *
bcm_CSRMatrixCreate( int num_rows,
		     int num_cols,
		     int num_nonzeros )
{
  bcm_CSRMatrix  *matrix;

  matrix = (bcm_CSRMatrix *) calloc(1, sizeof(bcm_CSRMatrix));

  bcm_CSRMatrixData(matrix) = NULL;
  bcm_CSRMatrixI(matrix)    = NULL;
  bcm_CSRMatrixJ(matrix)    = NULL;
  bcm_CSRMatrixNumRows(matrix) = num_rows;
  bcm_CSRMatrixNumCols(matrix) = num_cols;
  bcm_CSRMatrixNumNonzeros(matrix) = num_nonzeros;

  /* set defaults */
  bcm_CSRMatrixOwnsData(matrix) = 1;

  return matrix;
}
/*--------------------------------------------------------------------------
 * bcm_CSRMatrixDestroy
 *--------------------------------------------------------------------------*/

int 
bcm_CSRMatrixDestroy( bcm_CSRMatrix *matrix )
{
  int  ierr=0;

  if (matrix)
    {
      free(bcm_CSRMatrixI(matrix));
      if ( bcm_CSRMatrixOwnsData(matrix) )
	{
	  free(bcm_CSRMatrixData(matrix));
	  free(bcm_CSRMatrixJ(matrix));
	}
      free(matrix);
    }

  return ierr;
}

/*--------------------------------------------------------------------------
 * bcm_CSRMatrixInitialize
 *--------------------------------------------------------------------------*/

int 
bcm_CSRMatrixInitialize( bcm_CSRMatrix *matrix )
{
  int  num_rows     = bcm_CSRMatrixNumRows(matrix);
  int  num_nonzeros = bcm_CSRMatrixNumNonzeros(matrix);

  int  ierr=0;

  if ( ! bcm_CSRMatrixData(matrix) && num_nonzeros )
    bcm_CSRMatrixData(matrix) = (double *) calloc(num_nonzeros, sizeof(double));
  if ( ! bcm_CSRMatrixI(matrix) )
    bcm_CSRMatrixI(matrix)    = (int *) calloc(num_rows + 1, sizeof(int));
  if ( ! bcm_CSRMatrixJ(matrix) && num_nonzeros )
    bcm_CSRMatrixJ(matrix)    = (int *) calloc(num_nonzeros, sizeof(int));

  return ierr;
}

/*--------------------------------------------------------------------------
 * bcm_CSRMatrixSetDataOwner
 *--------------------------------------------------------------------------*/

int 
bcm_CSRMatrixSetDataOwner( bcm_CSRMatrix *matrix,
			   int              owns_data )
{
  int    ierr=0;

  bcm_CSRMatrixOwnsData(matrix) = owns_data;

  return ierr;
}

/*--------------------------------------------------------------------------
 * bcm_CSRMatrixRead
 *--------------------------------------------------------------------------*/

bcm_CSRMatrix *
bcm_CSRMatrixRead( char *file_name )
{
  bcm_CSRMatrix  *matrix;

  FILE    *fp;

  double  *matrix_data;
  int     *matrix_i;
  int     *matrix_j;
  int      num_rows;
  int      num_nonzeros;
  int      max_col = 0;

  int      file_base = 1;
   
  int      j;

  /*----------------------------------------------------------
   * Read in the data
   *----------------------------------------------------------*/

  fp = fopen(file_name, "r");

  fscanf(fp, "%d", &num_rows);

  matrix_i = (int *) calloc(num_rows + 1, sizeof(int));
  for (j = 0; j < num_rows+1; j++)
    {
      fscanf(fp, "%d", &matrix_i[j]);
      matrix_i[j] -= file_base;
    }

  num_nonzeros = matrix_i[num_rows];

  matrix = bcm_CSRMatrixCreate(num_rows, num_rows, matrix_i[num_rows]);
  bcm_CSRMatrixI(matrix) = matrix_i;
  bcm_CSRMatrixInitialize(matrix);

  matrix_j = bcm_CSRMatrixJ(matrix);
  for (j = 0; j < num_nonzeros; j++)
    {
      fscanf(fp, "%d", &matrix_j[j]);
      matrix_j[j] -= file_base;

      if (matrix_j[j] > max_col)
	{
	  max_col = matrix_j[j];
	}
    }

  matrix_data = bcm_CSRMatrixData(matrix);
  for (j = 0; j < matrix_i[num_rows]; j++)
    {
      fscanf(fp, "%le", &matrix_data[j]);
    }

  fclose(fp);

  bcm_CSRMatrixNumNonzeros(matrix) = num_nonzeros;
  bcm_CSRMatrixNumCols(matrix) = ++max_col;

  return matrix;
}

/*--------------------------------------------------------------------------
 * bcm_CSRMatrixPrint
 *--------------------------------------------------------------------------*/

int
bcm_CSRMatrixPrint( bcm_CSRMatrix *matrix,
		    char            *file_name )
{
  FILE    *fp;

  double  *matrix_data;
  int     *matrix_i;
  int     *matrix_j;
  int      num_rows;
  int      num_cols, nnz;
   
  int      file_base = 1;
   
  int      j;

  int      ierr = 0;

  /*----------------------------------------------------------
   * Print the matrix data
   *----------------------------------------------------------*/

  matrix_data = bcm_CSRMatrixData(matrix);
  matrix_i    = bcm_CSRMatrixI(matrix);
  matrix_j    = bcm_CSRMatrixJ(matrix);
  num_rows    = bcm_CSRMatrixNumRows(matrix);
  num_cols    = bcm_CSRMatrixNumCols(matrix);
  nnz         = bcm_CSRMatrixNumNonzeros(matrix);

  fp = fopen(file_name, "w");

  fprintf(fp, "%d  %d  \n", num_rows, nnz);

  for (j = 0; j <= num_rows; j++)
    {
      fprintf(fp, "%d\n", matrix_i[j] + file_base);
    }

  for (j = 0; j < matrix_i[num_rows]; j++)
    {
      fprintf(fp, "%d\n", matrix_j[j] + file_base);
    }

  if (matrix_data)
    {
      for (j = 0; j < matrix_i[num_rows]; j++)
	{
	  fprintf(fp, "%.14e\n", matrix_data[j]);
	}
    }
  else
    {
      fprintf(fp, "Warning: No matrix data!\n");
    }

  fclose(fp);

  return ierr;
}

/*--------------------------------------------------------------------------
 * bcm_CSRMatrixPrintMM : print a CSRMatrix in Matrix-Market format
 *--------------------------------------------------------------------------*/

int
bcm_CSRMatrixPrintMM( bcm_CSRMatrix *matrix,
		    char            *file_name )
{
  FILE    *fp;

  double  *matrix_data;
  int     *matrix_i;
  int     *matrix_j;
  int      num_rows;
  int      num_cols, nnz;
   
  int      file_base = 1;
   
  int      i,j;

  int      ierr = 0;

  /*----------------------------------------------------------
   * Print the matrix data
   *----------------------------------------------------------*/

  if (matrix == NULL) {
    fprintf(stderr,"Why do you want me to print a NULL matrix?\n");
    return(1);
  }
    
  matrix_data = bcm_CSRMatrixData(matrix);
  matrix_i    = bcm_CSRMatrixI(matrix);
  matrix_j    = bcm_CSRMatrixJ(matrix);
  num_rows    = bcm_CSRMatrixNumRows(matrix);
  num_cols    = bcm_CSRMatrixNumCols(matrix);
  nnz         = bcm_CSRMatrixNumNonzeros(matrix);

  fp = fopen(file_name, "w");
  fprintf(fp,"%s\n","%%MatrixMarket matrix coordinate real general");

  fprintf(fp, "%d  %d %d \n", num_rows, num_cols, nnz);

  for (i = 0; i < num_rows; i++)     {
    for (j=matrix_i[i]; j<matrix_i[i+1]; j++) {
      fprintf(fp, "%d   %d  %.18e\n", i+file_base,matrix_j[j] + file_base, matrix_data[j]);
    }
  }

  fclose(fp);

  return ierr;
}

/*--------------------------------------------------------------------------
 * bcm_CSRMatrixPrintHB: print a CSRMatrix in Harwell-Boeing format
 *--------------------------------------------------------------------------*/

int
bcm_CSRMatrixPrintHB( bcm_CSRMatrix *matrix_input,
		      char            *file_name )
{
  FILE            *fp;
  bcm_CSRMatrix *matrix;
  double          *matrix_data;
  int       *matrix_i;
  int       *matrix_j;
  int        num_rows;
  int        num_cols;
  int        file_base = 1;
  int        j, totcrd, ptrcrd, indcrd, valcrd, rhscrd;
  int        ierr = 0;

  /*----------------------------------------------------------
   * Print the matrix data
   *----------------------------------------------------------*/

  /* First transpose the input matrix, since HB is in CSC format */

  bcm_CSRMatrixTranspose(matrix_input, &matrix, 1);

  matrix_data = bcm_CSRMatrixData(matrix);
  matrix_i    = bcm_CSRMatrixI(matrix);
  matrix_j    = bcm_CSRMatrixJ(matrix);
  num_rows    = bcm_CSRMatrixNumRows(matrix);
  num_cols    = bcm_CSRMatrixNumCols(matrix);

  fp = fopen(file_name, "w");

  fprintf(fp, "%-70s  Key     \n", "Title");
  ptrcrd = num_rows;
  indcrd = matrix_i[num_rows];
  valcrd = matrix_i[num_rows];
  rhscrd = 0;
  totcrd = ptrcrd + indcrd + valcrd + rhscrd;
  fprintf (fp, "%14d%14d%14d%14d%14d\n",
	   totcrd, ptrcrd, indcrd, valcrd, rhscrd);
  fprintf (fp, "%-14s%14i%14i%14i%14i\n", "RUA",
	   num_rows, num_cols, valcrd, 0);
  fprintf (fp, "%-16s%-16s%-16s%26s\n", "(1I8)", "(1I8)", "(1E16.15)", "");

  for (j = 0; j <= num_rows; j++)
    {
      fprintf(fp, "%8d\n", matrix_i[j] + file_base);
    }

  for (j = 0; j < matrix_i[num_rows]; j++)
    {
      fprintf(fp, "%8d\n", matrix_j[j] + file_base);
    }

  if (matrix_data)
    {
      for (j = 0; j < matrix_i[num_rows]; j++)
	{
	  fprintf(fp, "%16.15e\n", matrix_data[j]);
	}
    }
  else
    {
      fprintf(fp, "Warning: No matrix data!\n");
    }

  fclose(fp);

  bcm_CSRMatrixDestroy(matrix);

  return ierr;
}

/*--------------------------------------------------------------------------
 * bcm_CSRMatrixCopy:
 * copys A to B, 
 * if copy_data = 0 only the structure of A is copied to B.
 * the routine does not check if the dimensions of A and B match !!! 
 *--------------------------------------------------------------------------*/

int 
bcm_CSRMatrixCopy( bcm_CSRMatrix *A, bcm_CSRMatrix *B, int copy_data )
{
  int  ierr=0;
  int  num_rows = bcm_CSRMatrixNumRows(A);
  int  num_cols = bcm_CSRMatrixNumCols(A);
  int  num_nonzeros = bcm_CSRMatrixNumNonzeros(A);
  int *A_i = bcm_CSRMatrixI(A);
  int *A_j = bcm_CSRMatrixJ(A);
  double *A_data;
  int *B_i = bcm_CSRMatrixI(B);
  int *B_j = bcm_CSRMatrixJ(B);
  double *B_data;

  int i, j;

  bcm_CSRMatrixNumRows(B)=num_rows;
  bcm_CSRMatrixNumCols(B)=num_cols;
  bcm_CSRMatrixNumNonzeros(B)=num_nonzeros;

  for (i=0; i < num_rows; i++)
    {
      B_i[i] = A_i[i];
      for (j=A_i[i]; j < A_i[i+1]; j++)
	{
	  B_j[j] = A_j[j];
	}
    }
  B_i[num_rows] = A_i[num_rows];
  if (copy_data)
    {
      A_data = bcm_CSRMatrixData(A);
      B_data = bcm_CSRMatrixData(B);
      for (i=0; i < num_rows; i++)
   	{
	  for (j=A_i[i]; j < A_i[i+1]; j++)
	    {
	      B_data[j] = A_data[j];
	    }
	}
    }
  return ierr;
}

/*--------------------------------------------------------------------------
 * bcm_CSRMatrixCloneStruct
 * Creates and returns a new copy of the argument, A.
 * Data is not copied, only structural information is reproduced.
 * Copying is a deep copy in that no pointers are copied; new arrays are
 * created where necessary.
 *--------------------------------------------------------------------------*/

bcm_CSRMatrix * bcm_CSRMatrixCloneStruct( bcm_CSRMatrix * A )
{
  int num_rows = bcm_CSRMatrixNumRows( A );
  int num_cols = bcm_CSRMatrixNumCols( A );
  int num_nonzeros = bcm_CSRMatrixNumNonzeros( A );
  bcm_CSRMatrix * B = bcm_CSRMatrixCreate( num_rows, num_cols, num_nonzeros );
  int * A_i;
  int * A_j;
  int * B_i;
  int * B_j;
  int i, j;

  bcm_CSRMatrixInitialize( B );

  A_i = bcm_CSRMatrixI(A);
  A_j = bcm_CSRMatrixJ(A);
  B_i = bcm_CSRMatrixI(B);
  B_j = bcm_CSRMatrixJ(B);

  for ( i=0; i<num_rows+1; ++i )  B_i[i] = A_i[i];
  for ( j=0; j<num_nonzeros; ++j )  B_j[j] = A_j[j];

  return B;
}

/*--------------------------------------------------------------------------
 * bcm_CSRMatrixClone
 * Creates and returns a new copy of the argument, A.
 * All data is copied.
 * Copying is a deep copy in that no pointers are copied; new arrays are
 * created where necessary.
 *--------------------------------------------------------------------------*/

bcm_CSRMatrix * bcm_CSRMatrixClone( bcm_CSRMatrix * A )
{
  if (A == NULL) {
    fprintf(stderr,"CSRMatrixClone: invalid A pointer\n" );
    return NULL;
  }
  int num_rows = bcm_CSRMatrixNumRows( A );
  int num_cols = bcm_CSRMatrixNumCols( A );
  int num_nonzeros = bcm_CSRMatrixNumNonzeros( A );
  int * A_i;
  int * A_j;
  double * A_data;
  int * B_i;
  int * B_j;
  double * B_data;
  int i, j;
  
  //fprintf(stderr,"Cloning: %d %d %d\n", num_rows, num_cols, num_nonzeros );
  bcm_CSRMatrix * B = bcm_CSRMatrixCreate( num_rows, num_cols, num_nonzeros );
  
  bcm_CSRMatrixInitialize( B );

  A_i = bcm_CSRMatrixI(A);
  A_j = bcm_CSRMatrixJ(A);
  A_data = bcm_CSRMatrixData(A);
  B_i = bcm_CSRMatrixI(B);
  B_j = bcm_CSRMatrixJ(B);
  B_data = bcm_CSRMatrixData(B);

  for ( i=0; i<num_rows+1; ++i )  B_i[i] = A_i[i];
  for ( j=0; j<num_nonzeros; ++j )  B_j[j] = A_j[j];
  for ( j=0; j<num_nonzeros; ++j )  B_data[j] = A_data[j];

  return B;
}
/*--------------------------------------------------------------------------
 * bcm_CSRMatrixTriU
 * Creates and returns a copy of the upper triangle of A starting from diagonal L.
 * L is 1,2 ...and represents the index of the upper diagonal (L=0 is the main diagonal)
 *--------------------------------------------------------------------------*/

bcm_CSRMatrix * bcm_CSRMatrixTriU( bcm_CSRMatrix * A, int L )
{
  int num_rows = bcm_CSRMatrixNumRows( A );
  int num_cols = bcm_CSRMatrixNumCols( A );
  int *A_i      = bcm_CSRMatrixI(A);
  int *A_j      = bcm_CSRMatrixJ(A);
  double *A_data   = bcm_CSRMatrixData(A);

  bcm_CSRMatrix *B;
  int * B_i;
  int * B_j;
  double *B_data;
  int j, k, i, l, num_nonzerosB, num_nnzrowB;

  num_nonzerosB=0;
  i=0;
  B_i = (int *) calloc(num_rows+1, sizeof(int));

  B_i[i]=0;
  for ( k=0; k< num_rows; ++k) 
    {
      num_nnzrowB=0;
      for (j=A_i[k]; j<A_i[k+1]; ++j)
	{
	  if(A_j[j] >= k+L) num_nnzrowB++;
	}

      num_nonzerosB+=num_nnzrowB;
      i++;
      B_i[i]=B_i[i-1]+num_nnzrowB;
    }
  B_i[num_rows]=num_nonzerosB;
   
  B = bcm_CSRMatrixCreate(num_rows, num_cols, num_nonzerosB);

  bcm_CSRMatrixInitialize(B);
  free(bcm_CSRMatrixI(B));
  bcm_CSRMatrixI(B) = B_i;
  B_j = bcm_CSRMatrixJ(B);
  B_data = bcm_CSRMatrixData(B);

  l=0;
  for (k=0; k<num_rows; ++k)
    {
    for(j=A_i[k]; j<A_i[k+1]; j++)
	{
	  if(A_j[j] >= k+L) 
	    {
	      B_j[l]=A_j[j];
	      B_data[l]=A_data[j];
	      l++;
	    }
	}

    }
  return B;
}
/*--------------------------------------------------------------------------
 * bcm_CSRMatrixTriL
 * Creates and returns a copy of the lower triangle of A starting from diagonal L.
 * L is 1,2 ...and it represents the index of the lower diagonal (L=0 is the main diagonal)
 *--------------------------------------------------------------------------*/

bcm_CSRMatrix *bcm_CSRMatrixTriL( bcm_CSRMatrix * A, int L )
{
  int num_rows = bcm_CSRMatrixNumRows( A );
  int num_cols = bcm_CSRMatrixNumCols( A );
  int *A_i      = bcm_CSRMatrixI(A);
  int *A_j      = bcm_CSRMatrixJ(A);
  double *A_data   = bcm_CSRMatrixData(A);

  bcm_CSRMatrix *B;
  int * B_i;
  int * B_j;
  double *B_data;
  int j, k, i, l, num_nonzerosB, num_nnzrowB;

  B_i = (int *) calloc(num_rows+1, sizeof(int));

  num_nonzerosB=0;
  i=0;
  B_i[i]=0;
  for ( k=0; k< num_rows; ++k) 
    {
      num_nnzrowB=0;
      for (j=A_i[k]; j<A_i[k+1]; ++j)
	{
	  if(A_j[j] <= k-L) num_nnzrowB=num_nnzrowB+1;
	}

      num_nonzerosB=num_nonzerosB+num_nnzrowB;
      i++;
      B_i[i]=B_i[i-1]+num_nnzrowB;
    }
  B_i[num_rows]=num_nonzerosB;

  B = bcm_CSRMatrixCreate(num_rows, num_cols, num_nonzerosB);

  bcm_CSRMatrixInitialize(B);
  free(bcm_CSRMatrixI(B));
  bcm_CSRMatrixI(B) = B_i;

  B_j = bcm_CSRMatrixJ(B);
  B_data = bcm_CSRMatrixData(B);

  l=0;
  for (k=0; k<num_rows; ++k)
    {
      for (j=A_i[k]; j<A_i[k+1]; ++j)
	{
	  if(A_j[j] <= k-L) 
	    {
	      B_j[l]=A_j[j];
	      B_data[l]=A_data[j];
	      l++;
	    }
	}
    }

  return B;
}
/*--------------------------------------------------------------------------
 * bcm_CSRMatrixDiag
 * Creates and returns a copy of the Diagonal of A
 *--------------------------------------------------------------------------*/

bcm_Vector * bcm_CSRMatrixDiag( bcm_CSRMatrix * A )
{
  int num_rows = bcm_CSRMatrixNumRows( A );
  int *A_i      = bcm_CSRMatrixI(A);
  int *A_j      = bcm_CSRMatrixJ(A);
  double *A_data   = bcm_CSRMatrixData(A);

  bcm_Vector *w;
  double *w_data;
  int j, k;

  w = bcm_VectorCreate(num_rows);
  bcm_VectorInitialize(w);

  w_data = bcm_VectorData(w);

  for (k=0; k<num_rows; ++k)
    {
      for (j=A_i[k]; j<A_i[k+1]; ++j)
	{
	  if(A_j[j] == k)  
	    {
	      w_data[k]=A_data[j];
	      if(w_data[k] < DBL_EPSILON) printf("Warning: NULL DIAGONAL ELEMENT IN DIAGONAL MATRIX\n");
	    }
	}
    }

  return w;
}
/*--------------------------------------------------------------------------
 * bcm_COO2CSRMatrixRead
 *--------------------------------------------------------------------------*/

bcm_CSRMatrix *
bcm_COO2CSRMatrixRead( char *file_name )
{
  bcm_CSRMatrix  *matrix;

  FILE    *fp;

  double  *matrix_value, *matrix_data;
  int     *matrix_cooi, *matrix_i;
  int     *matrix_cooj, *matrix_j;
  int      num_rows, num_cols;
  int      num_nonzeros;
  int      max_col = 0;

  int      file_base = 1;
   
  int      i, j, k, k0, iad;
  double   x;

  /*----------------------------------------------------------
   * Read in the data (matrix in COO format)
   *----------------------------------------------------------*/

  fp = fopen(file_name, "r");

  fscanf(fp, "%d", &num_rows);
  fscanf(fp, "%d", &num_nonzeros);

  matrix_cooi = (int *) calloc(num_nonzeros, sizeof(int));
  for (j = 0; j < num_nonzeros; j++)
    {
      fscanf(fp, "%d", &matrix_cooi[j]);
      matrix_cooi[j] -= file_base;
    }
  matrix_cooj = (int *) calloc(num_nonzeros, sizeof(int));
  for (j = 0; j < num_nonzeros; j++)
    {
      fscanf(fp, "%d", &matrix_cooj[j]);
      matrix_cooj[j] -= file_base;
      if (matrix_cooj[j] > max_col)
	{
	  max_col = matrix_cooj[j];
	}
    }
  matrix_value = (double *) calloc(num_nonzeros, sizeof(double));
  for (j = 0; j < num_nonzeros; j++) fscanf(fp, "%le", &matrix_value[j]);

  /*----------------------------------------------------------
   * Transform matrix from COO to CSR format
   *----------------------------------------------------------*/

  matrix_i = (int *) calloc(num_rows+1, sizeof(int));

  /* determine row lenght */
  for (j=0; j<num_nonzeros; j++) matrix_i[matrix_cooi[j]]=matrix_i[matrix_cooi[j]]+1;

  /* starting position of each row */
  k=0;
  for(j=0; j<= num_rows; j++)
    {
      k0=matrix_i[j];
      matrix_i[j]=k;
      k=k+k0;
    }
  matrix_j = (int *) calloc(num_nonzeros, sizeof(int));
  matrix_data = (double *) calloc(num_nonzeros, sizeof(double));

  /* go through the structure once more. Fill in output matrix */
  for(k=0; k<num_nonzeros; k++)
    {
      i=matrix_cooi[k];
      j=matrix_cooj[k];
      x=matrix_value[k];
      iad=matrix_i[i];
      matrix_data[iad]=x;
      matrix_j[iad]=j;
      matrix_i[i]=iad+1;
    }
  /* shift back matrix_i */
  for(j=num_rows-1; j>=0; j--) matrix_i[j+1]=matrix_i[j];
  matrix_i[0]=0;

  matrix = bcm_CSRMatrixCreate(num_rows, num_rows, num_nonzeros);
  bcm_CSRMatrixI(matrix) = matrix_i;
  bcm_CSRMatrixJ(matrix) = matrix_j;
  bcm_CSRMatrixData(matrix) = matrix_data;
  bcm_CSRMatrixNumNonzeros(matrix) = num_nonzeros;
  //bcm_CSRMatrixNumCols(matrix) = max_col++;

  free(matrix_cooi);
  free(matrix_cooj);
  free(matrix_value);
  fclose(fp);
  return matrix;
}
/*--------------------------------------------------------------------------
 * bcm_MM2CSRMatrixRead
 *--------------------------------------------------------------------------*/
#define BUFSIZE 1024
bcm_CSRMatrix *
bcm_MM2CSRMatrixRead( char *file_name )
{
  bcm_CSRMatrix  *matrix;

  FILE    *fp;
  char    banner[64], mtx[64], crd[64], data_type[64], storage_scheme[64];
  char    buffer[BUFSIZE+1];
  double  *matrix_value, *matrix_data, val;
  int     *matrix_cooi, *matrix_i;
  int     *matrix_cooj, *matrix_j;
  int      num_rows, num_cols, ri, cj;
  int      num_nonzeros, fr_nonzeros, allc_nonzeros;
  int      max_col = 0, is_general=0, is_symmetric=0;

  int      file_base = 1;
   
  int      i, j, k, k0, iad;
  double   x;

  /*----------------------------------------------------------
   * Read in the data (matrix in MM format)
   *----------------------------------------------------------*/

  fp = fopen(file_name, "r");

  fscanf(fp, "%s %s %s %s %s\n", banner, mtx, crd, data_type, storage_scheme);
  fgets(buffer,BUFSIZE,fp);
  for ( ; buffer[0]=='\%';  fgets(buffer,BUFSIZE,fp) );

  sscanf(buffer, "%d %d %d", &num_rows, &num_cols, &fr_nonzeros);

  if (strcmp(data_type,"real") !=0) {
    fprintf(stderr,"Error: we only read real matrices, not '%s'\n",data_type);
    fclose(fp);
    return(NULL);
  }
  
  if (strcmp(storage_scheme,"general")==0) {
    allc_nonzeros = fr_nonzeros;
    is_general=1;
  } else if (strcmp(storage_scheme,"symmetric")==0) {
    allc_nonzeros = 2*fr_nonzeros;
    is_symmetric=1;
  } else {
    fprintf(stderr,"Error: unhandled storage scheme '%s'\n",storage_scheme);
    fclose(fp);
    return(NULL);   
  }
  
  matrix_cooi = (int *) calloc(allc_nonzeros, sizeof(int));
  matrix_cooj = (int *) calloc(allc_nonzeros, sizeof(int));
  matrix_value = (double *) calloc(allc_nonzeros, sizeof(double));
  if (is_general) {
    num_nonzeros = fr_nonzeros;
    for (j = 0; j < fr_nonzeros; j++)
      {
	if (fgets(buffer,BUFSIZE,fp) != NULL) {
	  sscanf(buffer, "%d %d %le", &matrix_cooi[j], &matrix_cooj[j], &matrix_value[j]);
	  matrix_cooi[j] -= file_base;
	  matrix_cooj[j] -= file_base;	
	  if (matrix_cooj[j] > max_col)
	    {
	      max_col = matrix_cooj[j];
	    }
	} else {
	  fprintf(stderr,"Reading from MatrixMarket file failed\n");
	  fprintf(stderr,"Error while trying to read record %d of %d from file %s\n",
		  j,fr_nonzeros,file_name);
	  exit(-1);
	}
      }
  } else if (is_symmetric) {
    k = 0;
    for (j = 0; j < fr_nonzeros; j++)   {
      if (fgets(buffer,BUFSIZE,fp) != NULL) {
	sscanf(buffer, "%d %d %le", &ri, &cj, &val);
	ri -= file_base;
	cj -= file_base;
	if (cj > max_col)
	  max_col = cj;
	matrix_cooi[k]  = ri;
	matrix_cooj[k]  = cj;
	matrix_value[k] = val;
	k++;
	if (ri != cj) {
	  matrix_cooi[k]  = cj;
	  matrix_cooj[k]  = ri;
	  matrix_value[k] = val;
	  k++;
	}
      } else {
	fprintf(stderr,"Reading from MatrixMarket file failed\n");
	fprintf(stderr,"Error while trying to read record %d of %d from file %s\n",
		j,fr_nonzeros,file_name);
	fclose(fp);
	return(NULL);
      }
    }
    num_nonzeros = k;
  } else {
    fprintf(stderr,"Internal error: neither symmetric nor general ? \n");
    fclose(fp);
    return(NULL);
  }
  /*----------------------------------------------------------
   * Transform matrix from COO to CSR format
   *----------------------------------------------------------*/

  matrix_i = (int *) calloc(num_rows+1, sizeof(int));

  /* determine row lenght */
  for (j=0; j<num_nonzeros; j++) {
    if ((0<=matrix_cooi[j])&&(matrix_cooi[j]<num_rows)){
      matrix_i[matrix_cooi[j]]=matrix_i[matrix_cooi[j]]+1;
    } else {
      fprintf(stderr,"Wrong row index %d at position %d\n",matrix_cooi[j],j);
    }
  }

  /* starting position of each row */
  k=0;
  for(j=0; j<= num_rows; j++)
    {
      k0=matrix_i[j];
      matrix_i[j]=k;
      k=k+k0;
    }
  matrix_j = (int *) calloc(num_nonzeros, sizeof(int));
  matrix_data = (double *) calloc(num_nonzeros, sizeof(double));

  /* go through the structure once more. Fill in output matrix */
  for(k=0; k<num_nonzeros; k++)
    {
      i=matrix_cooi[k];
      j=matrix_cooj[k];
      x=matrix_value[k];
      iad=matrix_i[i];
      matrix_data[iad]=x;
      matrix_j[iad]=j;
      matrix_i[i]=iad+1;
    }
  /* shift back matrix_i */
  for(j=num_rows-1; j>=0; j--) matrix_i[j+1]=matrix_i[j];
  matrix_i[0]=0;

  matrix = bcm_CSRMatrixCreate(num_rows, num_cols, num_nonzeros);
  bcm_CSRMatrixI(matrix) = matrix_i;
  bcm_CSRMatrixJ(matrix) = matrix_j;
  bcm_CSRMatrixData(matrix) = matrix_data;
  bcm_CSRMatrixNumNonzeros(matrix) = num_nonzeros;
  //bcm_CSRMatrixNumCols(matrix) = max_col++;

  free(matrix_cooi);
  free(matrix_cooj);
  free(matrix_value);
  fclose(fp);
  return matrix;
}
/*--------------------------------------------------------------------------
 * bcm_CSRMatrixDiagScal
 * Apply diagonal scaling to matrix A: B=D^-1A 
 *--------------------------------------------------------------------------*/

bcm_CSRMatrix * bcm_CSRMatrixDiagScal( bcm_CSRMatrix * A, bcm_Vector *D )
{
  int num_rows = bcm_CSRMatrixNumRows( A );
  int num_cols = bcm_CSRMatrixNumCols( A );
  int num_nnzA = bcm_CSRMatrixNumNonzeros( A );
  int *A_i      = bcm_CSRMatrixI(A);
  int *A_j      = bcm_CSRMatrixJ(A);
  double *A_data   = bcm_CSRMatrixData(A);

  double *D_data   = bcm_CSRMatrixData(D);

  bcm_CSRMatrix *B;
  int * B_i;
  int * B_j;
  double *B_data;
  int j,k; 

  B = bcm_CSRMatrixCreate(num_rows, num_cols, num_nnzA);
  bcm_CSRMatrixInitialize(B);
  B_i = bcm_CSRMatrixI(B);
  B_j = bcm_CSRMatrixJ(B);
  B_data = bcm_CSRMatrixData(B);

  for ( k=0; k< num_rows; ++k) 
    {
      B_i[k]=A_i[k];
    }
  B_i[num_rows]=num_nnzA;
  for(k=0; k<num_nnzA; k++) 
   {
      B_data[k]=A_data[k];
      B_j[k]=A_j[k];
   }

  for (k=0; k<num_rows; ++k)
  {
      for (j=A_i[k]; j<A_i[k+1]; ++j) B_data[j]=A_data[j]/D_data[k];
  }

  return B;
}
