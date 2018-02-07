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

#include "bcm_amg.h"

/*----------------------------------------------------------------------------
 * bcm_CSRMatrixHMatch compute an half approximate maximum product matching 
 * in the graph of a sparse matrix by applying the algorithm in:
 * Preis, Linear Time 1/2-approximation algorithm for maximum weighted matching
 * in General Graph, in STACS'99, LNCS, vol.1563 (1999).
 *
 * Note: we expect a square matrix with symmetric sparsity pattern and null diagonal;
 *       Indeed we work with a strictly upper triangular matrix
 *--------------------------------------------------------------------------*/

int * bcm_CSRMatrixHMatch( bcm_CSRMatrix *B )
{

  int i, j, k, *rmatch;
  double  *c, alpha;
  int jbp, nzrows_B;
  double tmp=0.0;
  int rno, cno, nzrows_cno, startj, ljrowindex, ljjrowindex;
  int *jrowindex, *jjrowindex;

  bcm_CSRMatrix *BT; 

  int     *B_i = bcm_CSRMatrixI(B);
  int     *B_j = bcm_CSRMatrixJ(B);
  double  *B_data = bcm_CSRMatrixData(B);
  int  nrows_B  =  bcm_CSRMatrixNumRows(B);
  int  ncols_B  =  bcm_CSRMatrixNumCols(B);
  int  nnz_B  =  bcm_CSRMatrixNumNonzeros(B);

  assert(nrows_B==ncols_B);

  bcm_CSRMatrixTranspose ( B , &BT , 1); 

  int     *BT_i = bcm_CSRMatrixI(BT);
  int     *BT_j = bcm_CSRMatrixJ(BT);
  double  *BT_data = bcm_CSRMatrixData(BT);
  int  nrows_BT  =  bcm_CSRMatrixNumRows(BT); 

  nzrows_B=B_i[1]-B_i[0];

  /*prepare weight matrix for maximum product matching */

  /* I am working on B transpose since I need max along column of B 
-     so in CSR format is better working per row                                 */

  c = (double *) calloc(nrows_BT, sizeof(double));

  for(i=0; i<nrows_BT; i++) c[i]=-DBL_MIN;
  for(i=0; i<nrows_BT; ++i) 
    {
      for(j=BT_i[i]; j<BT_i[i+1]; ++j) 
	{
	  if(fabs(BT_data[j]) > 0.) 
	    {
	      tmp=log(fabs(BT_data[j]));
	      if(tmp >c[i]) c[i]=tmp;
	    }
	}
    }

  alpha=-DBL_MIN;
  for(i=0; i< nrows_BT; ++i)
    {
      for(j=BT_i[i]; j<BT_i[i+1]; j++)
	{
	  if(fabs(BT_data[j]) > 0.) 
	    {
	      tmp=c[i]-log(fabs(BT_data[j]));
	      if(tmp > alpha) alpha=tmp;
	    }
	}
    } 

  bcm_CSRMatrix * W= bcm_CSRMatrixCloneStruct(B);

  int     *W_i = bcm_CSRMatrixI(W);
  int     *W_j = bcm_CSRMatrixJ(W);
  int  nrows_W  =  bcm_CSRMatrixNumRows(W);
  double  *W_data=bcm_CSRMatrixData(W);

  /* needed for computing weights as in Bora code */

  double min_B=fabs(B_data[0]);
  for(i=1; i<nnz_B; i++)
    {
      if(min_B > fabs(B_data[i])) min_B=fabs(B_data[i]);
    } 

  /* end of weights as in Bora code */

  for(i=0; i< nrows_W; ++i)
    {
      for(j=W_i[i]; j<W_i[i+1]; ++j)
	{
	  /* As in Duff paper */
	            //if(fabs(B_data[j])>0.) W_data[j]=-log(c[B_j[j]])+log(fabs(B_data[j]));
		     // else W_data[j]=FLT_MAX; 
	  /* As in Scott paper for maximum cardinality */
	  // if(fabs(B_data[j])>0.) W_data[j]=alpha+log(fabs(B_data[j]))+(alpha-c[B_j[j]]);
	  //   else W_data[j]=FLT_MAX; 
	  /* As Bora code */
	  W_data[j]=log(fabs(B_data[j])/(0.999*min_B));  /* This seems to be generally better */
	}
    }

  rmatch = (int *) calloc(nrows_B, sizeof(int));

  for(i=0; i<nrows_B; ++i)  rmatch[i]=-1;

  jrowindex = (int *) calloc(nrows_B, sizeof(int));
  jjrowindex = (int *) calloc(nrows_B, sizeof(int));

  jbp=0;
  for(i=0; i< nrows_B; ++i)
    { 
      nzrows_B=B_i[i+1]-B_i[i];

      for(j=0; j<nzrows_B; ++j) jrowindex[j]=B_j[jbp+j];
      for(j=0; j<nzrows_B; ++j)
	{
	  rno=i;
	  cno=B_j[jbp+j];

	  startj=B_i[cno];
	  nzrows_cno=B_i[cno+1]-startj;

	  for(k=0; k<nzrows_cno; ++k) jjrowindex[k]=B_j[startj+k];

	  if(rmatch[rno] == -1 && rmatch[cno] == -1)
	    bcm_trymatch(rno,cno,W,jrowindex,nzrows_B,jjrowindex,nzrows_cno,rmatch);
        
	}

      jbp=jbp+nzrows_B;
    }
   
  bcm_CSRMatrixDestroy(BT);
  bcm_CSRMatrixDestroy(W);
  free(jrowindex);
  free(jjrowindex);
  free(c);
  
  return rmatch;
}


int bcm_trymatch(int rowindex, int colindex, bcm_CSRMatrix *W, int *jrowindex,
	     int ljrowindex, int *jcolindex, int ljcolindex, int *rmatch)
{

  int tryrowmatch, trycolmatch;
  int i, j, k, nzrow_W, startj, kindex;
  double cweight, nweight; 

  int     *W_i = bcm_CSRMatrixI(W);
  int     *W_j = bcm_CSRMatrixJ(W);
  int     nrows_W  =  bcm_CSRMatrixNumRows(W);
  double  *W_data=bcm_CSRMatrixData(W);

  k=-1;
  i=0;
  while (i<ljrowindex && k ==-1)
    {
      if(jrowindex[i]==colindex) k=i;
      i++;
    }
  if(k >= 0)
    {
      ljrowindex=ljrowindex-1;
      for(i=k; i<ljrowindex; ++i) jrowindex[i]=jrowindex[i+1];
    }

  k=-1;
  i=0;
  while(i<ljcolindex && k==-1)
    {
      if(jcolindex[i]==rowindex) k=i;
      i++;
    }
  if(k >= 0)
    {
      ljcolindex=ljcolindex-1;
      for(i=k; i<ljcolindex; ++i) jcolindex[i]=jcolindex[i+1];
    }

  nzrow_W=W_i[rowindex+1]-W_i[rowindex];
  startj=W_i[rowindex];
  kindex=-1;
  k=0;
  while(k<nzrow_W && kindex==-1) 
    {
      if(W_j[startj+k]==colindex) kindex=startj+k;
      k++;
    }
  cweight=W_data[kindex];

  while((rmatch[rowindex]==-1 && rmatch[colindex]==-1)
	&& (ljrowindex !=0 || ljcolindex !=0))    
    {

      if (rmatch[rowindex] == -1 && ljrowindex !=0)
	{
	  tryrowmatch=jrowindex[0];
	  nzrow_W=W_i[rowindex+1]-W_i[rowindex];
	  startj=W_i[rowindex];
	  k=0;
	  kindex=-1;
	  while(k<nzrow_W && kindex==-1) 
	    {
	      if(W_j[startj+k]==tryrowmatch) kindex=startj+k;
	      k++;
	    }
	  nweight=W_data[kindex];

	  ljrowindex=ljrowindex-1;
	  for(i=0; i<ljrowindex; ++i) jrowindex[i]=jrowindex[i+1];
        
	  if(nweight > cweight && rmatch[tryrowmatch]==-1)
	    {
          
              nzrow_W=W_i[tryrowmatch+1]-W_i[tryrowmatch];
	      int *trymatchindexrow;
	      trymatchindexrow= (int *) calloc(nzrow_W, sizeof(int));

	      startj=W_i[tryrowmatch];
	      for(k=0; k<nzrow_W; ++k) trymatchindexrow[k]=W_j[startj+k]; 

	      bcm_trymatch(rowindex,tryrowmatch,W,jrowindex,ljrowindex,
		       trymatchindexrow,nzrow_W,rmatch);
	      free(trymatchindexrow);
	    }
	}

      if(rmatch[colindex]==-1 && ljcolindex!=0)
	{
	  trycolmatch=jcolindex[0];
	  nzrow_W=W_i[colindex+1]-W_i[colindex];
	  startj=W_i[colindex];
	  k=0;
	  kindex=-1;
	  while(k<nzrow_W && kindex==-1) 
	    {
	      if(W_j[startj+k]==trycolmatch) kindex=startj+k;
	      k++;
	    }
	  nweight=W_data[kindex];

	  ljcolindex=ljcolindex-1;
	  for(i=0; i<ljcolindex; ++i) jcolindex[i]=jcolindex[i+1];
        
	  if(nweight > cweight && rmatch[trycolmatch]==-1)
	    {
	      nzrow_W=W_i[trycolmatch+1]-W_i[trycolmatch];
	      int *trymatchindexcol;
	      trymatchindexcol= (int *) calloc(nzrow_W, sizeof(int));

	      startj=W_i[trycolmatch];
	      for(k=0; k<nzrow_W; ++k) trymatchindexcol[k]=W_j[startj+k];

	      bcm_trymatch(colindex,trycolmatch,W,jcolindex,ljcolindex,
		       trymatchindexcol,nzrow_W,rmatch);
	      free(trymatchindexcol);
	    }
	}
    }
        
  if(rmatch[rowindex]==-1 & rmatch[colindex]==-1)
    {
      rmatch[rowindex]=colindex;
      rmatch[colindex]=rowindex;
    }

  return 0;
}
