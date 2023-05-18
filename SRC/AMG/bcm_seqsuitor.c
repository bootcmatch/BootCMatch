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
/* Suitor algorithm by Naim, Manne et al. 2015
 * contribution: Dario Pasquini, IAC-CNR March 2018 */
#include "bcm_amg.h"
#define true 1
#define false 0

int* suitor(int n, int nnz, int *row, int *col, double *val){

  int *vecsuitor, i, j;
  double *ws;
  
  vecsuitor = (int *) calloc(n, sizeof(int));
  ws = (double *) calloc(n, sizeof(double));

  for(i=0; i<n; i++){
    vecsuitor[i] = -1;
    ws[i] = -1;
  }
  
  // algorithm
  for(i=0; i<n; i++){
    int u = i;//row[i];
    int current = u;
    /*bool*/int done = false;

    while(!done){
      int partner = vecsuitor[current];
      double heaviest = ws[current];
      for(j=row[current]; j<row[current+1]; j++){
        int v = col[j];
        if(val[j] > heaviest && val[j] > ws[v]){
          partner = v;
          heaviest = val[j];
        }
      }
      done = true;

      if(heaviest != -1){
        int y = vecsuitor[partner];
        vecsuitor[partner] = current;
        ws[partner] = heaviest;
        if(y != -1){
          current = y;
          done = false;
        }
      }
    }
  }
  free(ws);
  return vecsuitor;
}
