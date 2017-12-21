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
#ifndef BCM_MATVEC_H_
#define BCM_MATVEC_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <ctype.h>

/******************************************************************************
 *
 * Basic data structures
 *
 *****************************************************************************/

/*--------------------------------------------------------------------------
 * CSR Matrix
 *--------------------------------------------------------------------------*/

typedef struct
{
   int     *i;
   int     *j;
   int      num_rows;
   int      num_cols;
   int      num_nonzeros;

   int      owns_data;

   double  *data;

} bcm_CSRMatrix;

/*--------------------------------------------------------------------------
 * Accessor functions for the CSR Matrix structure
 *--------------------------------------------------------------------------*/

#define bcm_CSRMatrixData(matrix)         ((matrix) -> data)
#define bcm_CSRMatrixI(matrix)            ((matrix) -> i)
#define bcm_CSRMatrixJ(matrix)            ((matrix) -> j)
#define bcm_CSRMatrixNumRows(matrix)      ((matrix) -> num_rows)
#define bcm_CSRMatrixNumCols(matrix)      ((matrix) -> num_cols)
#define bcm_CSRMatrixNumNonzeros(matrix)  ((matrix) -> num_nonzeros)
#define bcm_CSRMatrixOwnsData(matrix)     ((matrix) -> owns_data)


/******************************************************************************
 *
 * Vector
 *
 *****************************************************************************/

/*--------------------------------------------------------------------------
 * bcm_Vector
 *--------------------------------------------------------------------------*/

typedef struct
{
   double  *data;
   int      size;

   int      owns_data;

} bcm_Vector;

/*--------------------------------------------------------------------------
 * Accessor functions for the Vector structure
 *--------------------------------------------------------------------------*/

#define bcm_VectorData(vector)      ((vector) -> data)
#define bcm_VectorSize(vector)      ((vector) -> size)
#define bcm_VectorOwnsData(vector)  ((vector) -> owns_data)

/* csr_matop.c  */
bcm_CSRMatrix *bcm_CSRMatrixAdd ( bcm_CSRMatrix *A , bcm_CSRMatrix *B );
bcm_CSRMatrix *bcm_CSRMatrixMultiply ( bcm_CSRMatrix *A , bcm_CSRMatrix *B );
bcm_CSRMatrix *bcm_CSRMatrixDeleteZeros ( bcm_CSRMatrix *A , double tol );
int bcm_CSRMatrixTranspose ( bcm_CSRMatrix *A , bcm_CSRMatrix **AT , int data );
int bcm_CSRMatrixReorder ( bcm_CSRMatrix *A );
double bcm_CSRMatrixSumElts ( bcm_CSRMatrix *A );
double bcm_CSRMatrixInfNorm ( bcm_CSRMatrix *A );
double bcm_VectorANorm( bcm_CSRMatrix *A, bcm_Vector *vector );
bcm_Vector * bcm_CSRMatrixLSolve(bcm_CSRMatrix *L, bcm_Vector *D, bcm_Vector *f);
bcm_Vector * bcm_CSRMatrixUSolve(bcm_CSRMatrix *U, bcm_Vector *D, bcm_Vector *f);
int bcm_CSRMatrixSort ( bcm_CSRMatrix *A );
int bcm_SortMatrixRow(int nz, int *ja, double *val);


/* csr_matrix.c */
bcm_CSRMatrix *bcm_CSRMatrixCreate ( int num_rows , int num_cols , int num_nonzeros );
int bcm_CSRMatrixDestroy ( bcm_CSRMatrix *matrix );
int bcm_CSRMatrixInitialize ( bcm_CSRMatrix *matrix );
int bcm_CSRMatrixSetDataOwner ( bcm_CSRMatrix *matrix , int owns_data );
bcm_CSRMatrix *bcm_CSRMatrixRead ( char *file_name );
bcm_CSRMatrix *bcm_COO2CSRMatrixRead ( char *file_name );
bcm_CSRMatrix *bcm_MM2CSRMatrixRead( char *file_name );
int bcm_CSRMatrixPrint ( bcm_CSRMatrix *matrix , char *file_name );
int bcm_CSRMatrixPrintHB ( bcm_CSRMatrix *matrix_input , char *file_name );
int bcm_CSRMatrixPrintMM( bcm_CSRMatrix *matrix, char            *file_name );
int bcm_CSRMatrixCopy ( bcm_CSRMatrix *A , bcm_CSRMatrix *B , int copy_data );
bcm_CSRMatrix * bcm_CSRMatrixCloneStruct ( bcm_CSRMatrix *A );
bcm_CSRMatrix * bcm_CSRMatrixClone ( bcm_CSRMatrix *A );
bcm_CSRMatrix * bcm_CSRMatrixTriU( bcm_CSRMatrix * A, int L );
bcm_CSRMatrix * bcm_CSRMatrixTriL( bcm_CSRMatrix * A, int L );
bcm_Vector * bcm_CSRMatrixDiag( bcm_CSRMatrix * A);
bcm_CSRMatrix * bcm_CSRMatrixDiagScal( bcm_CSRMatrix * A, bcm_Vector *D );


/* csr_matvec.c */
int bcm_CSRMatrixMatvec ( double alpha , bcm_CSRMatrix *A , bcm_Vector *x , double beta , bcm_Vector *y );
int bcm_CSRMatrixMatvecT ( double alpha , bcm_CSRMatrix *A , bcm_Vector *x , double beta , bcm_Vector *y );

/* vector.c */
bcm_Vector *bcm_VectorCreate ( int size );
int bcm_VectorFreeData ( bcm_Vector *vector );
int bcm_VectorDestroy ( bcm_Vector *vector );
int bcm_VectorInitialize ( bcm_Vector *vector );
bcm_Vector *bcm_VectorRead ( char *file_name );
bcm_Vector *bcm_VectorMMRead ( char *file_name );
int bcm_VectorPrint ( bcm_Vector *vector , char *file_name );
int bcm_VectorSetConstantValues ( bcm_Vector *v , double value );
int bcm_VectorSetRandomValues ( bcm_Vector *v , int seed );
int bcm_VectorCopy ( bcm_Vector *x , bcm_Vector *y );
bcm_Vector *bcm_VectorCloneDeep ( bcm_Vector *x );
bcm_Vector *bcm_VectorCloneShallow ( bcm_Vector *x );
int bcm_VectorScale ( double alpha , bcm_Vector *y );
int bcm_VectorAxpy ( double alpha , bcm_Vector *x , bcm_Vector *y );
double bcm_VectorInnerProd ( bcm_Vector *x , bcm_Vector *y );
double bcm_VectorSumElts ( bcm_Vector *vector );
double bcm_VectorNorm( bcm_Vector *vector );
bcm_Vector * bcm_VectorDiagScal( bcm_Vector *u, bcm_Vector * D );

/* csr_relax.c */
int  bcm_CSRMatrixRelax( bcm_CSRMatrix *A, bcm_CSRMatrix *L, bcm_CSRMatrix *U, bcm_Vector *D, bcm_Vector *f, int relax_type, double relax_weight, bcm_Vector *u); 

#endif
