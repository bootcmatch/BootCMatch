/*
                BootCMatch
     Bootstrap AMG based on Compatible Matching version 1.0
    (C) Copyright 2017
                       Pasqua D'Ambra    IAC-CNR
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

/******************************************************************************
 *
 * functions for bcm_Vector structure
 *
 *****************************************************************************/

#include "bcm_matvec.h"

/*--------------------------------------------------------------------------
 * bcm_VectorCreate
 *--------------------------------------------------------------------------*/

bcm_Vector *
bcm_VectorCreate( int size )
{
  bcm_Vector  *vector;

  vector = (bcm_Vector *) calloc(1, sizeof(bcm_Vector));

  bcm_VectorData(vector) = NULL;
  bcm_VectorSize(vector) = size;

  /* set defaults */
  bcm_VectorOwnsData(vector) = 1;

  return vector;
}

/*--------------------------------------------------------------------------
 * bcm_VectorDestroy
 *--------------------------------------------------------------------------*/

int
bcm_VectorDestroy( bcm_Vector *vector )
{
   int  ierr=0;

   if (vector)
   {
      if ( bcm_VectorOwnsData(vector) )
      {
         free(bcm_VectorData(vector));
      }
      free(vector);
   }

   return ierr;
}

/*--------------------------------------------------------------------------
 * bcm_VectorFreeData
 *--------------------------------------------------------------------------*/

int 
bcm_VectorFreeData( bcm_Vector *vector )
{
  int  ierr=0;

  if (vector)
    {
      if ( bcm_VectorOwnsData(vector) )
	{
	  free(bcm_VectorData(vector));
	}
    }

  return ierr;
}

/*--------------------------------------------------------------------------
 * bcm_VectorInitialize
 *--------------------------------------------------------------------------*/

int 
bcm_VectorInitialize( bcm_Vector *vector )
{
  int  size = bcm_VectorSize(vector);
  int  ierr = 0;

  if ( ! bcm_VectorData(vector) )
    bcm_VectorData(vector) = (double *) calloc(size, sizeof(double));

    ++ierr;

  return ierr;
}

/*--------------------------------------------------------------------------
 * bcm_VectorSetDataOwner
 *--------------------------------------------------------------------------*/

int 
bcm_VectorSetDataOwner( bcm_Vector *vector,
			int           owns_data   )
{
  int    ierr=0;

  bcm_VectorOwnsData(vector) = owns_data;

  return ierr;
}

/*--------------------------------------------------------------------------
 * bcm_VectorRead
 *--------------------------------------------------------------------------*/

bcm_Vector *
bcm_VectorRead( char *file_name )
{
  bcm_Vector  *vector;

  FILE    *fp;

  double  *data;
  int      size;
   
  int      j;

  /*----------------------------------------------------------
   * Read in the data
   *----------------------------------------------------------*/

  fp = fopen(file_name, "r");

  fscanf(fp, "%d", &size);

  vector = bcm_VectorCreate(size);
  bcm_VectorInitialize(vector);

  data = bcm_VectorData(vector);
  for (j = 0; j < size; j++)
    {
      fscanf(fp, "%le", &data[j]);
    }

  fclose(fp);

  return vector;
}
/*--------------------------------------------------------------------------
 * bcm_VectorMMRead: Read Vector from Matrix Market file
 *--------------------------------------------------------------------------*/

bcm_Vector *
bcm_VectorMMRead( char *file_name )
{
  bcm_Vector  *vector;

  FILE    *fp;
  char    banner[64], mtx[64], crd[64], data_type[64], storage_scheme[64];

  double  *data;
  int      size1, size2, idx;
   
  int      j;

  /*----------------------------------------------------------
   * Read in the data
   *----------------------------------------------------------*/

  fp = fopen(file_name, "r");

  fscanf(fp, "%s %s %s %s %s\n", banner, mtx, crd, data_type, storage_scheme);
  fscanf(fp, "%d %d", &size1, &size2);

  vector = bcm_VectorCreate(size1);
  bcm_VectorInitialize(vector);

  data = bcm_VectorData(vector);
  for (j = 0; j < size1; j++)
    {
      fscanf(fp, "%le\n", &data[j]);
    }

  fclose(fp);

  return vector;
}

/*--------------------------------------------------------------------------
 * bcm_VectorPrint
 *--------------------------------------------------------------------------*/

int
bcm_VectorPrint( bcm_Vector *vector,
		 char         *file_name )
{
  FILE    *fp;

  int ierr=0;

  fp = fopen(file_name, "w");
  bcm_VectorPrintfp(vector,fp)  ;
  fclose(fp);

  return ierr;
}


int
bcm_VectorPrintfp( bcm_Vector *vector,
		   FILE    *fp)
{

  double  *data;
  int      size;
   
  int      i;

  int      ierr = 0;

  /*----------------------------------------------------------
   * Print in the data
   *----------------------------------------------------------*/

  data = bcm_VectorData(vector);
  size = bcm_VectorSize(vector);


      fprintf(fp, "%d\n", size);

      for (i = 0; i < size; i++)
	{
	  fprintf(fp, "%.14e\n", data[i]);
	}


  return ierr;
}

/*--------------------------------------------------------------------------
 * bcm_VectorSetConstantValues
 *--------------------------------------------------------------------------*/

int
bcm_VectorSetConstantValues( bcm_Vector *v,
			     double        value )
{
  double  *vector_data = bcm_VectorData(v);
  int      size        = bcm_VectorSize(v);
           
  int      i;
           
  int      ierr  = 0;

  for (i = 0; i < size; i++)
    vector_data[i] = value;

  return ierr;
}

/*--------------------------------------------------------------------------
 * bcm_VectorSetRandomValues
 *
 *     returns vector of values randomly distributed between -1.0 and +1.0
 *--------------------------------------------------------------------------*/

int
bcm_VectorSetRandomValues( bcm_Vector *v,
			   int           seed )
{
  double  *vector_data = bcm_VectorData(v);
  int      size        = bcm_VectorSize(v);

  void bcm_SeedRand (int seed);
  double bcm_Rand();
           
  int      i;
           
  int      ierr  = 0;
  bcm_SeedRand(seed);

  for (i = 0; i < size; i++)
    vector_data[i] = 2.0 * bcm_Rand() - 1.0;

  return ierr;
}

/*--------------------------------------------------------------------------
 * bcm_VectorCopy
 * copies data from x to y
 * y should have already been initialized at the same size as x
 *--------------------------------------------------------------------------*/

int
bcm_VectorCopy( bcm_Vector *x,
		bcm_Vector *y )
{
  double  *x_data = bcm_VectorData(x);
  double  *y_data = bcm_VectorData(y);
  int      size   = bcm_VectorSize(x);
           
  int      i;
           
  int      ierr = 0;

  for (i = 0; i < size; i++)
    y_data[i] = x_data[i];

  return ierr;
}

/*--------------------------------------------------------------------------
 * bcm_VectorCloneDeep
 * Returns a complete copy of x - a deep copy, with its own copy of the data.
 *--------------------------------------------------------------------------*/

bcm_Vector *
bcm_VectorCloneDeep( bcm_Vector *x )
{
  int      size   = bcm_VectorSize(x);
  bcm_Vector * y = bcm_VectorCreate( size);

  bcm_VectorInitialize(y);
  bcm_VectorCopy( x, y );

  return y;
}

/*--------------------------------------------------------------------------
 * bcm_VectorCloneShallow
 * Returns a complete copy of x - a shallow copy, pointing the data of x
 *--------------------------------------------------------------------------*/

bcm_Vector *
bcm_VectorCloneShallow( bcm_Vector *x )
{
  int      size   = bcm_VectorSize(x);
  bcm_Vector * y = bcm_VectorCreate( size );

  bcm_VectorData(y) = bcm_VectorData(x);
  bcm_VectorSetDataOwner( y, 0 );
  bcm_VectorInitialize(y);

  return y;
}

/*--------------------------------------------------------------------------
 * bcm_VectorScale
 *--------------------------------------------------------------------------*/

int
bcm_VectorScale( double        alpha,
		 bcm_Vector *y     )
{
  double  *y_data = bcm_VectorData(y);
  int      size   = bcm_VectorSize(y);
           
  int      i;
           
  int      ierr = 0;

  for (i = 0; i < size; i++)
    y_data[i] *= alpha;

  return ierr;
}

/*--------------------------------------------------------------------------
 * bcm_VectorAxpy
 *--------------------------------------------------------------------------*/

int
bcm_VectorAxpy( double        alpha,
		bcm_Vector *x,
		bcm_Vector *y     )
{
  double  *x_data = bcm_VectorData(x);
  double  *y_data = bcm_VectorData(y);
  int      size   = bcm_VectorSize(x);
           
  int      i;
           
  int      ierr = 0;

  for (i = 0; i < size; i++)
    y_data[i] += alpha * x_data[i];

  return ierr;
}

/*--------------------------------------------------------------------------
 * bcm_VectorInnerProd
 *--------------------------------------------------------------------------*/

double   bcm_VectorInnerProd( bcm_Vector *x,
			      bcm_Vector *y )
{
  double  *x_data = bcm_VectorData(x);
  double  *y_data = bcm_VectorData(y);
  int      size   = bcm_VectorSize(x);
           
  int      i;

  double      result = 0.0;

  for (i = 0; i < size; i++)
    result += y_data[i] * x_data[i];

  return result;
}

/*--------------------------------------------------------------------------
 * bcm_VectorSumElts:
 * Returns the sum of all vector elements.
 *--------------------------------------------------------------------------*/

double bcm_VectorSumElts( bcm_Vector *vector )
{
  double sum = 0;
  double * data = bcm_VectorData( vector );
  int size = bcm_VectorSize( vector );
  int i;

  for ( i=0; i<size; ++i ) sum += data[i];

  return sum;
}
/*--------------------------------------------------------------------------
 * bcm_VectorNorm:
 * Returns the Euclidean norm of the vector.
 *--------------------------------------------------------------------------*/

double bcm_VectorNorm( bcm_Vector *vector )
{
  double norm = 0;

  norm=bcm_VectorInnerProd(vector,vector);
  norm=sqrt(norm);

  return norm;
}

/*--------------------------------------------------------------------------
 * bcm_VectorDiagScal
 * Apply diagonal scaling to vector v: w=D^-1v 
 *--------------------------------------------------------------------------*/

bcm_Vector * bcm_VectorDiagScal(bcm_Vector *u, bcm_Vector * D )
{
  int num_rows = bcm_VectorSize( D );
  double *D_data   = bcm_VectorData(D);
  double *u_data   = bcm_VectorData(u);

  bcm_Vector *w;
  double *w_data;
  int k; 
  w = bcm_VectorCreate(num_rows);
  bcm_VectorInitialize(w);
  w_data = bcm_VectorData(w);

  for (k=0; k<num_rows; ++k)
    {
      w_data[k]=u_data[k]/D_data[k];
    }

  return w;
}
