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
/******************************************************************************
 *
 * functions for CSR sparse matrix-vector operations
 *
 *****************************************************************************/

#include "bcm_matvec.h"

/*--------------------------------------------------------------------------
 * bcm_CSRMatrixMatvec: axpy operation
 *
 *   Performs y <- alpha * A * x + beta * y
 *--------------------------------------------------------------------------*/

int
bcm_CSRMatrixMatvec( double           alpha,
		     bcm_CSRMatrix *A,
		     bcm_Vector    *x,
		     double           beta,
		     bcm_Vector    *y     )
{
  double     *A_data   = bcm_CSRMatrixData(A);
  int        *A_i      = bcm_CSRMatrixI(A);
  int        *A_j      = bcm_CSRMatrixJ(A);
  int         num_rows = bcm_CSRMatrixNumRows(A);
  int         num_cols = bcm_CSRMatrixNumCols(A);

  double     *x_data = bcm_VectorData(x);
  double     *y_data = bcm_VectorData(y);
  int         x_size = bcm_VectorSize(x);
  int         y_size = bcm_VectorSize(y);

  double      temp, tempx;

  int         i, j, jj;

  int         m;

  int        ierr = 0;
 
  double dlamch_(char *cmach);


  if (num_cols != x_size)
    ierr = 1;

  if (num_rows != y_size)
    ierr = 2;

  if (num_cols != x_size && num_rows != y_size)
    ierr = 3;

  if (fabs(alpha) <= DBL_EPSILON)
    {
      if (fabs(beta) <= DBL_EPSILON)
	{
	  for (i = 0; i < num_rows; i++)
	    y_data[i] =0.0;
	}
      else
	{
	  for (i = 0; i < num_rows; i++)
	    y_data[i] *= beta;
	}     
      
      return ierr;
    }

  /*-----------------------------------------------------------------------
   * y = (beta/alpha)*y
   *-----------------------------------------------------------------------*/

  if (fabs(beta) < DBL_EPSILON)
    beta = 0.0;

  /*-----------------------------------------------------------------
   * y += A*x
   *-----------------------------------------------------------------*/

	  if (alpha == 1.0)
	    {
	      if (beta == 0.0)
		{
		  for (i = 0; i < num_rows; i++)
		    {
		      temp = 0.0;
		      for (jj = A_i[i]; jj < A_i[i+1]; jj++)
			temp += A_data[jj] * x_data[A_j[jj]];
		      y_data[i] = temp;
		    }
		}
	      else if (beta == 1.0)
		{
		  for (i = 0; i < num_rows; i++)
		    {
		      temp = y_data[i];
		      for (jj = A_i[i]; jj < A_i[i+1]; jj++)
			temp += A_data[jj] * x_data[A_j[jj]];
		      y_data[i] = temp;
		    }
		}
	      else if (beta == -1.0)
		{
		  for (i = 0; i < num_rows; i++)
		    {
		      temp =  - y_data[i];
		      for (jj = A_i[i]; jj < A_i[i+1]; jj++)
			temp += A_data[jj] * x_data[A_j[jj]];
		      y_data[i] = temp;
		    }
		}
	      else
		{
		  for (i = 0; i < num_rows; i++)
		    {
		      temp = 0.0;
		      for (jj = A_i[i]; jj < A_i[i+1]; jj++)
			temp += A_data[jj] * x_data[A_j[jj]];
		      y_data[i] = temp+beta*y_data[i] ;
		    }
		}
	    }
	  else if (alpha == -1.0)
	    {
	      if (beta == 0.0)
		{
		  for (i = 0; i < num_rows; i++)
		    {
		      temp = 0.0;
		      for (jj = A_i[i]; jj < A_i[i+1]; jj++)
			temp -= A_data[jj] * x_data[A_j[jj]];
		      y_data[i] = temp;
		    }
		}
	      else if (beta == 1.0)
		{
		  for (i = 0; i < num_rows; i++)
		    {
		      temp =  y_data[i];
		      for (jj = A_i[i]; jj < A_i[i+1]; jj++)
			temp -= A_data[jj] * x_data[A_j[jj]];
		      y_data[i] = temp;
		    }
		}
	      else if (beta == -1.0)
		{
		  for (i = 0; i < num_rows; i++)
		    {
		      temp =  - y_data[i];;
		      for (jj = A_i[i]; jj < A_i[i+1]; jj++)
			temp -= A_data[jj] * x_data[A_j[jj]];
		      y_data[i] = temp ;
		    }
		}
	      else
		{
		  for (i = 0; i < num_rows; i++)
		    {
		      temp = 0.0;
		      for (jj = A_i[i]; jj < A_i[i+1]; jj++)
			temp -= A_data[jj] * x_data[A_j[jj]];
		      y_data[i] = temp+beta*y_data[i] ;
		    }
		}
	    }
	  else 
	    {
	      if (beta == 0.0)
		{
		  for (i = 0; i < num_rows; i++)
		    {
		      temp = 0.0;
		      for (jj = A_i[i]; jj < A_i[i+1]; jj++)
			temp += A_data[jj] * x_data[A_j[jj]];
		      y_data[i] = alpha*temp;
		    }
		}
	      else
		{
		  for (i = 0; i < num_rows; i++)
		    {
		      temp = 0.0;
		      for (jj = A_i[i]; jj < A_i[i+1]; jj++)
			temp += A_data[jj] * x_data[A_j[jj]];
		      y_data[i] = alpha*temp+beta*y_data[i] ;
		    }
		}
	    }

  return ierr;
}

/*--------------------------------------------------------------------------
 * bcm_CSRMatrixMatvecT: axpy operation involving sparse matrix transpose
 *
 *   Performs y <- alpha * A^T * x + beta * y
 *
 *--------------------------------------------------------------------------*/

int
bcm_CSRMatrixMatvecT( double           alpha,
		      bcm_CSRMatrix *A,
		      bcm_Vector    *x,
		      double           beta,
		      bcm_Vector    *y     )
{
  double     *A_data    = bcm_CSRMatrixData(A);
  int        *A_i       = bcm_CSRMatrixI(A);
  int        *A_j       = bcm_CSRMatrixJ(A);
  int         num_rows  = bcm_CSRMatrixNumRows(A);
  int         num_cols  = bcm_CSRMatrixNumCols(A);

  double     *x_data = bcm_VectorData(x);
  double     *y_data = bcm_VectorData(y);
  int         x_size = bcm_VectorSize(x);
  int         y_size = bcm_VectorSize(y);
  double      temp;

  int         i, i1, j, jv, jj, ns, ne, size, rest;

  int         ierr  = 0;

  double dlamch_(char *cmach);

  if (num_rows != x_size)
    ierr = 1;

  if (num_cols != y_size)
    ierr = 2;

  if (num_rows != x_size && num_cols != y_size)
    ierr = 3;

  if (fabs(alpha) <= DBL_EPSILON)
    {
      for (i = 0; i < num_cols; i++)
	y_data[i] *= beta;

      return ierr;
    }

  /*-----------------------------------------------------------------------
   * y = (beta/alpha)*y
   *-----------------------------------------------------------------------*/

  temp = beta / alpha;
   
  if (temp != 1.0)
    {
      if (fabs(temp) <= DBL_EPSILON)
	{
	  for (i = 0; i < num_cols; i++)
	    y_data[i] = 0.0;
	}
      else
	{
	  for (i = 0; i < num_cols; i++)
	    y_data[i] *= temp;
	}
    }

  /*-----------------------------------------------------------------
   * y += A^T*x
   *-----------------------------------------------------------------*/
  for (i = 0; i < num_rows; i++)
    {
	  for (jj = A_i[i]; jj < A_i[i+1]; jj++)
            {
	      j = A_j[jj];
	      y_data[j] += A_data[jj] * x_data[i];
            }
    }
  /*-----------------------------------------------------------------
   * y = alpha*y
   *-----------------------------------------------------------------*/

  if (alpha != 1.0)
    {
      for (i = 0; i < num_cols; i++)
	y_data[i] *= alpha;
    }

  return ierr;
}
