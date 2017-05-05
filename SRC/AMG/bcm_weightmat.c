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

#include "bcm_matvec.h"

/*----------------------------------------------------------------------------
 * bcm_CSRMatrixAhat computes the matrix Ahat to which apply the weighted matching 
 * starting from A and w. 
 * Matrix A is considered s.p.d. (symmetric with nonzero diagonal entries is sufficient)
 * However, we can also modify it for working on unsymmetric matrices.
 *--------------------------------------------------------------------------*/

bcm_CSRMatrix *
bcm_CSRMatrixAhat( bcm_CSRMatrix *A,
		   bcm_Vector *w )
{

  int     *A_i = bcm_CSRMatrixI(A);
  int     *A_j = bcm_CSRMatrixJ(A);
  double  *A_data = bcm_CSRMatrixData(A);

  int     sizew    =  bcm_VectorSize(w);
  double  *w_data  =  bcm_VectorData(w);

  bcm_CSRMatrix  *B, *AH;
  double         *AH_data;
  int            *AH_i;
  int            *AH_j;

  int           k, l, icol, irow, nnz_AH, i, j;
  int           L=1;
  double        Aedge11, Aedge12, Aedge22;
  double        wedge0, wedge1, normwedge;
/* for old weights */
//  double        wedgeperp0, wedgeperp1, wtmp1, wtmp2, Aedge21;

  double dlamch_(char *cmach);

  // TriU is the correct one, since Preis matching is implemented working on upper triangle. 
  B=bcm_CSRMatrixTriU(A, L );

  int  nrows_B  =  bcm_CSRMatrixNumRows(B);
  int  ncols_B  =  bcm_CSRMatrixNumCols(B);
  int  nnz_B  =  bcm_CSRMatrixNumNonzeros(B);
  int  *B_i = bcm_CSRMatrixI(B);
  int  *B_j = bcm_CSRMatrixJ(B);

  assert(nrows_B==sizew);

  AH= bcm_CSRMatrixCreate(nrows_B, ncols_B, nnz_B);
  bcm_CSRMatrixInitialize(AH);

  AH_i = bcm_CSRMatrixI(AH);
  AH_j = bcm_CSRMatrixJ(AH);
  AH_data= bcm_CSRMatrixData(AH);

  nnz_AH=0;
  for(k=0; k< nrows_B; ++k)
    {
      irow=k;

      AH_i[k]=B_i[k];

      for(l=B_i[k]; l< B_i[k+1]; l++)
      {
          icol=B_j[l];

	  AH_j[l]=icol;

           /* new weights coming from D-orthogonality condition */

	  if(irow != icol)
	  {
	      wedge0=w_data[irow];
	      wedge1=w_data[icol];

             for(i=A_i[irow]; i<A_i[irow+1]; i++)
             {
               if(A_j[i]==irow) Aedge11=A_data[i];
               else if(A_j[i]==icol) Aedge12=A_data[i];
             }
             for(i=A_i[icol]; i<A_i[icol+1]; i++)
             {
               if(A_j[i]==icol) Aedge22=A_data[i];
             }
	     normwedge=Aedge11*pow(wedge0,2)+Aedge22*pow(wedge1,2);
	     if(normwedge > DBL_EPSILON)
             {
	        AH_data[nnz_AH]=1.0-((2.0*wedge0*wedge1*Aedge12)/normwedge);
	        nnz_AH++;
             }
	     else
	     {
	       AH_data[nnz_AH]=DBL_EPSILON;
	       nnz_AH++;
             }
	  } 
/* end of computation of new weights */

/* computation of old weights (standard L2 orthogonality condition) */
   /*    if(irow != icol)
       {
         wedge0=w_data[irow];
         wedge1=w_data[icol];

         normwedge=sqrt(pow(wedge0,2)+pow(wedge1,2));
       
         if(normwedge > DBL_EPSILON)
         {
          wedgeperp0=-w_data[icol]/normwedge;
          wedgeperp1=w_data[irow]/normwedge;

          for(i=A_i[irow]; i<A_i[irow+1]; i++)
          {
            if(A_j[i]==irow) Aedge11=A_data[i];
            else if(A_j[i]==icol) Aedge21=A_data[i];
          }

          for(i=A_i[icol]; i<A_i[icol+1]; i++)
          {
            if(A_j[i]==irow) Aedge12=A_data[i];
            else if(A_j[i]==icol) Aedge22=A_data[i];
          }


          wtmp1=wedgeperp0*Aedge11+wedgeperp1*Aedge21;
          wtmp2=wedgeperp0*Aedge12+wedgeperp1*Aedge22;

          AH_data[nnz_AH]=wtmp1*wedgeperp0+wtmp2*wedgeperp1;
          nnz_AH++;
        }
        else
        {
          AH_data[nnz_AH]=DBL_EPSILON;
          nnz_AH++;
        }
      } */
/* end of computation of old weights */
      }
      Aedge11=0.0;
      Aedge22=0.0;
      Aedge12=0.0;
      wedge0=0;
      wedge1=0;
    } 
    AH_i[nrows_B]=B_i[nrows_B];

    assert(nnz_AH==nnz_B);
   
    bcm_CSRMatrixDestroy(B);

    return AH;
}
