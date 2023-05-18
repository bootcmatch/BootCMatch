/*
                BootCMatch
     Bootstrap AMG based on Compatible Matching version 1.0
    (C) Copyright 2017
                       Pasqua D'Ambra    IAC-CNR
                       Panayot S. Vassilevski Portland State University, OR, USA

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

/* ***********************************************************************************
 *
 *  functions for building AMG hierarchies using Compatible
 *  Weighted (Max Product) Matching
 *
 *********************************************************************************** */

#include "bcm_amg.h"
#ifdef HAVE_HSL
#include "hsl_mc64d.h"
#endif
#ifdef HAVE_SUPERLU
#include "slu_ddefs.h"
#endif
#ifdef HAVE_SPRAL
#include "spral.h"
#endif
/* ***********************************************************************************
 *
 *  bcm_CSRMatchingAgg computes prolongator matrix by using Matching Vector
 *
 *********************************************************************************** */
#define MATCH_HSL      1
#define MATCH_SPRAL    2
#define MATCH_SUITOR   3
#define MATCH_PREIS    0
bcm_CSRMatrix * bcm_CSRMatchingAgg(bcm_CSRMatrix *A, bcm_Vector **w,
				   bcm_CSRMatrix **P, int match_type, int num_sweeps,
				   int max_sizecoarse, int max_levels, int *ftcoarse,
				   int cr_it, int cr_relax_type, double cr_relax_weight)
{ 
  bcm_CSRMatrix *A_temp, At, **A_tmp, **P_tmp, *P_temp , **L_tmp, **U_tmp, *R, *Ac, *PMM, *Pagg;
  bcm_Vector **D_tmp;
  bcm_Vector  *w_temp, *w_temp1, *rhs;
  int lev, i, sizecoarse, k, real_num_sweeps, nsize_w, nsize_A;
  double coarseratio;
  double timematching;
  double timegalerkin;

  bcm_AMGHierarchy *amg_hierarchy_tmp; /* aux Hierarchy */

  nsize_w=bcm_VectorSize(*w);
  /* initialize rhs vector for relaxations on homegeneous systems */
  rhs=bcm_VectorCreate(nsize_w);
  bcm_VectorInitialize(rhs);

  /* initialize wtemp for restriction of the current smooth vector at each level */
  w_temp=bcm_VectorCreate(nsize_w);
  bcm_VectorInitialize(w_temp);
  w_temp1=bcm_VectorCreate(nsize_w);
  bcm_VectorInitialize(w_temp1);
  bcm_VectorCopy(*w,w_temp);
  nsize_w=bcm_VectorSize(*w); 
  nsize_A=bcm_CSRMatrixNumRows(A); 
  assert(nsize_A == nsize_w);

  amg_hierarchy_tmp= bcm_AMGHierarchyCreate(num_sweeps+1);
  bcm_AMGHierarchyInitialize(amg_hierarchy_tmp);

  lev=0;
  sizecoarse=nsize_A;
  A_tmp= bcm_AMGHierarchyAArray(amg_hierarchy_tmp); 
  P_tmp= bcm_AMGHierarchyPArray(amg_hierarchy_tmp);
  L_tmp= bcm_AMGHierarchyLArray(amg_hierarchy_tmp);
  U_tmp= bcm_AMGHierarchyUArray(amg_hierarchy_tmp);
  D_tmp= bcm_AMGHierarchyDArray(amg_hierarchy_tmp);

  At=*A;
  A_tmp[lev]=&At;
  L_tmp[lev]=bcm_CSRMatrixTriL(A,1);
  U_tmp[lev]=bcm_CSRMatrixTriU(A,1);
  D_tmp[lev]=bcm_CSRMatrixDiag(A);
  int ierr;

  fprintf(stderr,"MatchingPairAgg num_sweeps: %d\n",num_sweeps);
  /* cycle for composition of pairwise prolongator to obtain more aggressive coarsening */
  double timematchingtot=0.0;
  double timegalerkintot=0.0;
  for(i=1; i<=num_sweeps; i++){
    /* build prolongator by pairwise aggregation based on compatible weighted matching */
    timematching=time_getWallclockSeconds(); 
    ierr=bcm_CSRMatchingPairAgg(A_tmp[i-1],w_temp,&Pagg,match_type);
    timematching=time_getWallclockSeconds()-timematching;
    timematchingtot=timematchingtot+timematching;
    fprintf(stderr,"From MatchingPairAgg: %d %p\n",ierr,Pagg);
    if (ierr==0) {
	P_tmp[i-1] = Pagg;
	bcm_CSRMatrixTranspose(Pagg, &R, 1);

	/* build coarse matrix */
        timegalerkin=time_getWallclockSeconds(); 
	A_temp=bcm_CSRMatrixMultiply(R,A_tmp[i-1]);
	A_tmp[i]=bcm_CSRMatrixMultiply(A_temp,P_tmp[i-1]);
	sizecoarse=bcm_CSRMatrixNumRows(A_tmp[i]);
        timegalerkin=time_getWallclockSeconds()-timegalerkin;
        timegalerkin=timegalerkintot+timegalerkin;
	bcm_CSRMatrixDestroy(A_temp);
	/* relax restricted w  */
	bcm_VectorSize(w_temp1)=sizecoarse;
	bcm_CSRMatrixMatvec(1.0,R,w_temp,0.0,w_temp1);
	bcm_VectorSize(rhs)=sizecoarse;
	/* extract triangular and diagonal part of coarse matrices */
	L_tmp[i]=bcm_CSRMatrixTriL(A_tmp[i],1);
	U_tmp[i]=bcm_CSRMatrixTriU(A_tmp[i],1);
	D_tmp[i]=bcm_CSRMatrixDiag(A_tmp[i]);
	bcm_CSRMatrixDestroy(R);
	for(k=1; k <= cr_it; ++k) 
	  bcm_CSRMatrixRelax(A_tmp[i], L_tmp[i],
			     U_tmp[i], D_tmp[i],
			     rhs, cr_relax_type, cr_relax_weight, w_temp1); 
	bcm_VectorSize(w_temp)=sizecoarse; 
	bcm_VectorCopy(w_temp1,w_temp);
	bcm_VectorSize(*w)=sizecoarse; 
	bcm_VectorCopy(w_temp,*w);
	sizecoarse  = bcm_CSRMatrixNumRows(A_tmp[i]);
	coarseratio = bcm_CSRMatrixNumRows(A_tmp[i-1]);
	coarseratio = coarseratio/sizecoarse;

      }else {
      break;
    }
    /* If sizecoarse meets required criteria exit the loop */
    if (coarseratio <= 1.2) *ftcoarse=100;
    if (sizecoarse<=(*ftcoarse)*max_sizecoarse) {
      i++;
      break;
    } 
  }
  real_num_sweeps=i - 1;
  printf("Time for aggregation sweeps:  %e\n", timematching);
  printf("Time for Galerkin product:  %e\n", timegalerkin);

  PMM = bcm_CSRMatrixClone(P_tmp[0]);
  if (real_num_sweeps > 1) {
    for (i=1; i< real_num_sweeps; i++) {
      P_temp = PMM;
      PMM=bcm_CSRMatrixMultiply(P_temp,P_tmp[i]);
      bcm_CSRMatrixDestroy(P_temp);
    }
  }
  *P= PMM;

  Ac=bcm_CSRMatrixClone(A_tmp[real_num_sweeps]);

  A_tmp[0]=NULL;
  for (i=0; i< real_num_sweeps; i++){
    bcm_CSRMatrixDestroy(A_tmp[i]);
    bcm_CSRMatrixDestroy(L_tmp[i]);
    bcm_CSRMatrixDestroy(U_tmp[i]);
    bcm_CSRMatrixDestroy(P_tmp[i]);
    bcm_VectorDestroy(D_tmp[i]); 
    A_tmp[i]=NULL;
    L_tmp[i]=NULL;
    U_tmp[i]=NULL;
    P_tmp[i]=NULL;
    D_tmp[i]=NULL;
  }
  bcm_AMGHierarchyDestroy(amg_hierarchy_tmp);
  bcm_VectorDestroy(rhs); 
  bcm_VectorDestroy(w_temp); 
  bcm_VectorDestroy(w_temp1); 
  return Ac;
}

int 
bcm_CSRMatchingPairAgg(bcm_CSRMatrix *A, bcm_Vector *w, bcm_CSRMatrix **P, int match_type)

{ 
  int sizew = bcm_VectorSize(w);
  double *w_data = bcm_VectorData(w);

  int ierr=0, i, j, k, nzeros;
  int ncolc=0, *markc, npairs=0, nsingle=0;
  double *wtempc, wagg0, wagg1, normwagg;

  bcm_CSRMatrix *AHT;
  bcm_CSRMatrix *AH_full;
  int *ia, *jcp;
  double *val;

#if defined(HAVE_HSL) || defined(HAVE_SPRAL)
  int nrows_L;
  double *rscaling, *cscaling;
#endif
#if defined(HAVE_HSL) 
  struct mc64_control mc_cntrl;
  struct mc64_info    mc_info;  
#endif
#if  defined(HAVE_SPRAL)
  struct spral_scaling_auction_options options;
  struct spral_scaling_auction_inform inform;
#endif

  int *p;
  int nrows_A = bcm_CSRMatrixNumRows(A);
  assert(nrows_A == sizew);


  /* build matrix Ahat from A */

  bcm_CSRMatrix *AH = bcm_CSRMatrixAhat(A,w,match_type);

  /* build permutation vector by weighted matching  */


  switch (match_type) {
#ifdef HAVE_HSL
  /* call to MC64 for optimal matching */
  case MATCH_HSL:
    /* Here we are using the fact that the upper triangle in CSR
       is exactly identical in memory to the lower triangle in CSC */
    nzeros = bcm_CSRMatrixNumNonzeros(AH);
    nrows_L = bcm_CSRMatrixNumRows(AH);
    bcm_CSRMatrixTranspose(AH, &AHT, 1);
    AH_full=bcm_CSRMatrixAdd(AH,AHT);

    jcp = bcm_CSRMatrixI(AH_full);
    ia  = bcm_CSRMatrixJ(AH_full);
    val = bcm_CSRMatrixData(AH_full);

    /* Now go for MC64 */
    mc64_default_control(&mc_cntrl);
    p   = (int *)calloc(2*nrows_L,sizeof(int));
    int job   = 5;  /*maximize the product of the weight matrix */
    int mtype = 2;
    fprintf(stderr,"Calling mc64\n");
    mc64_matching(job,mtype,nrows_L,nrows_L,jcp,ia,val,&mc_cntrl,&mc_info,p,NULL);

    for(i=0; i<nrows_L; ++i) p[i]=p[nrows_L+i];

    if (mc_info.flag < 0) {
      fprintf(stderr,"Warning: problem from mc64_matching %d\n",mc_info.flag);
      exit(1);
    }
    else if (mc_info.flag == 1)
      fprintf(stderr,"Warning: matrix appears structurally singular  %d\n",mc_info.flag);
    bcm_CSRMatrixDestroy(AHT);
    bcm_CSRMatrixDestroy(AH_full);
    break;
#endif
#if defined(HAVE_SPRAL)
  /* call to SPRAL for auction-based matching */
  case MATCH_SPRAL:
    nzeros = bcm_CSRMatrixNumNonzeros(AH);
    nrows_L = bcm_CSRMatrixNumRows(AH);
    bcm_CSRMatrixTranspose(AH, &AHT, 1);
    AH_full=bcm_CSRMatrixAdd(AH,AHT);
    
    jcp = bcm_CSRMatrixI(AH_full);
    ia  = bcm_CSRMatrixJ(AH_full);
    val = bcm_CSRMatrixData(AH_full);
    
    
    p   = (int *) calloc(2*nrows_L,sizeof(int));
    rscaling = (double *)calloc(nrows_L,sizeof(double));
    cscaling = (double *)calloc(nrows_L,sizeof(double));
    fprintf(stderr,"Calling spral_scaling_auction_unsym\n");
    spral_scaling_auction_default_options(&options);
    spral_scaling_auction_unsym(nrows_L, nrows_L, jcp, ia, val,
				rscaling, cscaling, p,
				&options, &inform);
    if(inform.flag<0) {
      printf("spral_scaling_auction_unsym() returned with error %5d", inform.flag);
      exit(1);
    }
    
    for(i=0; i<nrows_L; ++i) {
      p[nrows_L+i] = p[i];
    }
    bcm_CSRMatrixDestroy(AHT);
    bcm_CSRMatrixDestroy(AH_full);
    free(rscaling);
    free(cscaling);      
    break;
#endif
  /* call function for half-approximate matching based on Preis algorithm*/
  case MATCH_PREIS:
    fprintf(stderr,"Calling PREIS matching\n");
    p = bcm_CSRMatrixHMatch(AH);
    break;
  /* call function for half-approximate matching based on Suitor algorithm*/
  case MATCH_SUITOR:
  /*prepare weight matrix for maximum product matching */
    bcm_CSRMatrixTranspose(AH, &AHT, 1);
    AH_full=bcm_CSRMatrixAdd(AH,AHT);
    nzeros = bcm_CSRMatrixNumNonzeros(AH_full);

    jcp = bcm_CSRMatrixI(AH_full);
    ia  = bcm_CSRMatrixJ(AH_full);
    val = bcm_CSRMatrixData(AH_full);
    p=suitor(nrows_A, nzeros, jcp, ia, val);
    bcm_CSRMatrixDestroy(AHT);
    bcm_CSRMatrixDestroy(AH_full);
    break;
  default:
    fprintf(stderr,"Error: unknown matching algorithm\n");
    return(-1);
  }

  markc = (int *) calloc(nrows_A, sizeof(int));
  for(i=0; i<nrows_A; ++i) markc[i]=-1;

  wtempc = (double *) calloc(nrows_A, sizeof(double));


  for(i=0; i<nrows_A; ++i)
    {
      j=p[i];

      if((j>=0) && (i != j))
	{
	  if(markc[i] == -1 && markc[j]== -1)
	    {
	      wagg0=w_data[i];
	      wagg1=w_data[j];
	      normwagg=sqrt(pow(wagg0,2)+pow(wagg1,2));
//	      using the De-norm
	      //normwagg=sqrt(D_data[i]*pow(wagg0,2)+D_data[j]*pow(wagg1,2));
	      if(normwagg > DBL_EPSILON)
		{
		  markc[i]=ncolc;
		  markc[j]=ncolc;
		  wtempc[i]=w_data[i]/normwagg;
		  wtempc[j]=w_data[j]/normwagg;
		  ncolc++;
		  npairs++;
		}
	    }
	}
    }

  for(i=0; i<nrows_A; ++i)
    {
      if(markc[i]==-1)
	{
	  if(fabs(w_data[i]) <= DBL_EPSILON)
	    {
	      /* only-fine grid node: corresponding null row in the prolongator */
	      markc[i]=ncolc-1;
	      wtempc[i]=0.0;
	    }
	  else
	    {
	      markc[i]=ncolc;
	      ncolc++;
	      wtempc[i]=w_data[i]/fabs(w_data[i]);
	      nsingle++;
	    }
	}
    }

  int ncoarse=npairs+nsingle;
  fprintf(stderr,"npairs: %d nsingle: %d\n",npairs,nsingle);
  

  assert(ncolc == ncoarse);

  /* Each row of P has only 1 element. It can be zero in case
     of only-fine grid variable or singleton */

  *P=bcm_CSRMatrixCreate(nrows_A,ncolc,nrows_A);
  bcm_CSRMatrixInitialize(*P);

  if (ncoarse > 0)
    {
      int *P_i=bcm_CSRMatrixI(*P);
      int *P_j=bcm_CSRMatrixJ(*P);
      double *P_data=bcm_CSRMatrixData(*P);

      for(i=0; i<nrows_A; i++)
	{
	  P_i[i]=i;
	  P_j[i]=markc[i];
	  P_data[i]=wtempc[i];
	}
      P_i[nrows_A]=nrows_A;

     }

  bcm_CSRMatrixDestroy(AH);
  free(p);
  free(markc);
  free(wtempc);
  return 0;
}

/* ***********************************************************************************
 *
 *  bcm_AdaptiveCoarsening computes hierarchy of operators starting
 *  from a given smooth vector and using aggregation based on compatible weighted matching
 *
 *********************************************************************************** */

bcm_AMGHierarchy * bcm_AdaptiveCoarsening(bcm_AMGBuildData *amg_data)
{

  bcm_CSRMatrix *A, *P, *A_temp, *R, *P_temp,**P2;
  bcm_Vector *w, *wc, *w_temp, *w_temp1;
  int nsize_A, nsize_w;
  double cmplxfinal, cmplxsmoothed;
  double wcmplxfinal, wcmplxsmoothed;
  bcm_Vector *rhs;
  int cr_relax_type, cr_it,  coarse_solver;
  double cr_relax_weight;
  int i, lev, sizecoarse, max_levels, max_sizecoarse, j, k;
  double dlamch_(char *cmach);
  double  normw, normold, normnew, ratioAc;
  double coarseratio=2.0, avcoarseratio=0.0;

  bcm_CSRMatrix **A_array, **P_array, **L_array, **U_array ;
  bcm_Vector **D_array;

  bcm_AMGHierarchy *amg_hierarchy; /* Hierarchy to be built */
  A=bcm_AMGBuildDataCSRMatrix(amg_data);
  w=bcm_AMGBuildDataSmoothVector(amg_data);

  max_levels=bcm_AMGBuildDataMaxLevels(amg_data);

  nsize_A=bcm_CSRMatrixNumRows(A); 
  nsize_w=bcm_VectorSize(w); 


  assert(nsize_A == nsize_w);

  cr_relax_type   = bcm_AMGBuildDataCRRelaxType(amg_data);
  cr_relax_weight = bcm_AMGBuildDataCRRelaxWeight(amg_data);
  cr_it           = bcm_AMGBuildDataCRIterations(amg_data);
  max_levels      = bcm_AMGBuildDataMaxLevels(amg_data);
  max_sizecoarse  = bcm_AMGBuildDataMaxCoarseSize(amg_data);
  coarse_solver   = bcm_AMGBuildDataCoarseSolver(amg_data);

  /* initialize rhs vector for relaxations on homegeneous systems */
  rhs=bcm_VectorCreate(nsize_w);
  bcm_VectorInitialize(rhs);

  /* initialize wtemp for restriction of the current smooth vector at each level */
  w_temp=bcm_VectorCreate(nsize_w);
  bcm_VectorInitialize(w_temp);
  w_temp1=bcm_VectorCreate(nsize_w);
  bcm_VectorInitialize(w_temp1);

  bcm_VectorCopy(w,w_temp);
  int num_sweeps=bcm_AMGBuildDataSweepNumber(amg_data)+1;
  amg_hierarchy= bcm_AMGHierarchyCreate(max_levels);
  bcm_AMGHierarchyInitialize(amg_hierarchy);
  normw=bcm_VectorANorm(A,w_temp);
  lev=0;
  sizecoarse=nsize_A;

  A_array= bcm_AMGHierarchyAArray(amg_hierarchy); 
  P_array= bcm_AMGHierarchyPArray(amg_hierarchy);
  L_array= bcm_AMGHierarchyLArray(amg_hierarchy);
  U_array= bcm_AMGHierarchyUArray(amg_hierarchy);
  D_array= bcm_AMGHierarchyDArray(amg_hierarchy);

  double timesmoother;
  double timesmoothertot=0.0;

  timesmoother=time_getWallclockSeconds();
  A_array[lev]=A;
  A_temp=A_array[lev];
  L_array[lev]=bcm_CSRMatrixTriL(A,1);
  U_array[lev]=bcm_CSRMatrixTriU(A,1);
  D_array[lev]=bcm_CSRMatrixDiag(A);
  timesmoother=time_getWallclockSeconds()-timesmoother;
  fprintf(stderr,"Time for smoother setup: %e\n",timesmoother);
  timesmoothertot=timesmoothertot+timesmoother;


  /* Because of this initialization, and since MatrixDestroy
     & friends correctly handle NULL pointers, we are allowed to
     call XXXDestroy prior to any call that returns a new object;
     hence we can safeguard against memory leaks */
  for (i=lev+1; i<max_levels; i++) {
    A_array[i]=NULL;
    P_array[i]=NULL;
    L_array[i]=NULL;
    U_array[i]=NULL;
    D_array[i]=NULL;
  }
  double *tmp_do;
  /* relax on w_temp */
  for(i=1; i <= cr_it; ++i)
    bcm_CSRMatrixRelax(A_array[0], L_array[0], U_array[0], D_array[0], rhs,
		       cr_relax_type, cr_relax_weight, w_temp);

  cmplxsmoothed=bcm_CSRMatrixNumNonzeros(A);
  wcmplxsmoothed=bcm_CSRMatrixNumNonzeros(A);

  int match_type=bcm_AMGBuildDataAggMatchType(amg_data);

  double time1=time_getWallclockSeconds();
  int ftcoarse=1;
  int real_num_sweeps;
  int ierr;

  avcoarseratio=0.0;
  if (normw > DBL_EPSILON)   {
    j = 1;
    while(lev < max_levels && sizecoarse>ftcoarse*max_sizecoarse)  {

/* truncation of the smooth vector */
/*  if(lev == 0) {
  double maxw=0.0;
  double *w_data;
  w_data=bcm_VectorData(w_temp);
  for(i=0; i<nsize_A; ++i)
  {
   if(fabs(w_data[i])>maxw) maxw=fabs(w_data[i]);  
  }
  for(i=0; i<nsize_A; ++i) w_data[i]=w_data[i]/maxw;  
  for(i=0; i<nsize_A; ++i) 
  {
     if(fabs(w_data[i]) < 0.001) 
     {
      w_data[i]=0.0;  
     }
  } 
  } */
/* end of truncation process */
      
      A_array[j]=bcm_CSRMatchingAgg(A_array[j-1],&w_temp,&P,match_type, 
				    num_sweeps, max_sizecoarse, max_levels, &ftcoarse,
				    cr_it, cr_relax_type, cr_relax_weight);
      
      P_array[j-1]=P;
      sizecoarse=bcm_CSRMatrixNumRows(A_array[j]);
      coarseratio=bcm_CSRMatrixNumRows(A_array[j-1]);
      coarseratio=coarseratio/sizecoarse;
      avcoarseratio=avcoarseratio+coarseratio;
  timesmoother=time_getWallclockSeconds();
      L_array[j]=bcm_CSRMatrixTriL(A_array[j],1);
      U_array[j]=bcm_CSRMatrixTriU(A_array[j],1);
      D_array[j]=bcm_CSRMatrixDiag(A_array[j]);
  timesmoother=time_getWallclockSeconds()-timesmoother;
  timesmoothertot=timesmoothertot+timesmoother;
  fprintf(stderr,"Time for smoother setup: %e\n",timesmoother);
      j++;  
      lev=j-1;
    }
    lev=lev+1;
    bcm_AMGHierarchyNumLevels(amg_hierarchy)=lev;
  } else {
    printf("Warning: no need to build multigrid since the matrix is well conditioned\n");
  }

  cmplxfinal=bcm_CSRMatrixNumNonzeros(A);
  for(i=1; i< lev; i++){
    cmplxfinal += bcm_CSRMatrixNumNonzeros(A_array[i]); 
  }  
  cmplxfinal=cmplxfinal/bcm_CSRMatrixNumNonzeros(A);
  bcm_AMGHierarchyOpCmplx(amg_hierarchy)=cmplxfinal;

  wcmplxfinal=bcm_CSRMatrixNumNonzeros(A);
  for(i=1; i< lev; i++){
    wcmplxfinal += pow(2,i)*bcm_CSRMatrixNumNonzeros(A_array[i]); 
  }  
  wcmplxfinal=wcmplxfinal/bcm_CSRMatrixNumNonzeros(A);
  bcm_AMGHierarchyOpCmplxW(amg_hierarchy)=wcmplxfinal;

  avcoarseratio=avcoarseratio/(lev-1);

  bcm_AMGHierarchyAvgCratio(amg_hierarchy)=avcoarseratio;

  double time2=time_getWallclockSeconds()-time1;

  printf("Time for coarsening:  %e\n", time2);
  printf("Time for smoother setup: %e\n",timesmoothertot);

  /* free unused memory */
  for(i=lev; i<max_levels; i++) bcm_CSRMatrixDestroy(P_array[i]);
  for(i=lev; i<=max_levels; i++) bcm_CSRMatrixDestroy(A_array[i]);
  for(i=lev; i<=max_levels; i++) bcm_CSRMatrixDestroy(U_array[i]);
  for(i=lev; i<=max_levels; i++) bcm_CSRMatrixDestroy(L_array[i]);
  for(i=lev; i<=max_levels; i++) bcm_VectorDestroy(D_array[i]);

  bcm_VectorDestroy(rhs); 
  bcm_VectorDestroy(w_temp); 
  bcm_VectorDestroy(w_temp1); 

  /* apply one sweep of a weighted-Jacobi smoother to the prolongators
   * for smoothed-type aggregation */

  time1=time_getWallclockSeconds();
  if(bcm_AMGBuildDataAggInterpType(amg_data))
    {
      int nsize, nnzatemp;
      double *Atemp_data, omega=0.0;
      int *Atemp_i, *Atemp_j;
      for(j=0; j<lev-1; j++)
	{
	  nsize=bcm_CSRMatrixNumRows(A_array[j]);

	  A_temp=bcm_CSRMatrixDiagScal(A_array[j],D_array[j]);
	  omega=4.0/(3.0*bcm_CSRMatrixInfNorm(A_temp));

	  nnzatemp=bcm_CSRMatrixNumNonzeros(A_temp);
	  Atemp_data=bcm_CSRMatrixData(A_temp);
	  Atemp_i=bcm_CSRMatrixI(A_temp);
	  Atemp_j=bcm_CSRMatrixJ(A_temp);
	  for(k=0; k<nnzatemp; k++) Atemp_data[k]=-omega*Atemp_data[k];

	  for(i=0; i<nsize; i++)
            {
	      for(k=Atemp_i[i]; k<Atemp_i[i+1]; k++) 
		if(Atemp_j[k]==i) Atemp_data[k]=1.0+Atemp_data[k];
            }

          P_temp=bcm_CSRMatrixCreate(bcm_CSRMatrixNumRows(P_array[j]),
				     bcm_CSRMatrixNumCols(P_array[j]),
				     bcm_CSRMatrixNumNonzeros(P_array[j]));
          P_temp = bcm_CSRMatrixClone(P_array[j]);
	  if (P_array[j]!= NULL) bcm_CSRMatrixDestroy(P_array[j]);
	  P_array[j]=bcm_CSRMatrixMultiply(A_temp,P_temp);
	  bcm_CSRMatrixDestroy(P_temp);
	  bcm_CSRMatrixTranspose(P_array[j],&R,1);

	  if (A_temp!= NULL) bcm_CSRMatrixDestroy(A_temp);
	  A_temp=bcm_CSRMatrixMultiply(R,A_array[j]);

	  //if (A_array[j+1]!= NULL) bcm_CSRMatrixDestroy(A_array[j+1]);
	  A_array[j+1]=bcm_CSRMatrixMultiply(A_temp,P_array[j]);
	  if (L_array[j+1]!= NULL) bcm_CSRMatrixDestroy(L_array[j+1]);
	  if (U_array[j+1]!= NULL) bcm_CSRMatrixDestroy(U_array[j+1]);
	  if (D_array[j+1]!= NULL) bcm_VectorDestroy(D_array[j+1]);
	  L_array[j+1]=bcm_CSRMatrixTriL(A_array[j+1],1);
	  U_array[j+1]=bcm_CSRMatrixTriU(A_array[j+1],1);
	  D_array[j+1]=bcm_CSRMatrixDiag(A_array[j+1]);

	  bcm_CSRMatrixDestroy(R); 
	  bcm_CSRMatrixDestroy(A_temp);
	}
         cmplxsmoothed=bcm_CSRMatrixNumNonzeros(A);
         for(i=1; i< lev; i++){
            cmplxsmoothed += bcm_CSRMatrixNumNonzeros(A_array[i]); 
          }  
         cmplxsmoothed=cmplxsmoothed/bcm_CSRMatrixNumNonzeros(A);
         bcm_AMGHierarchyOpCmplx(amg_hierarchy)=cmplxsmoothed;

         wcmplxfinal=bcm_CSRMatrixNumNonzeros(A);
         for(i=1; i< lev; i++){
           wcmplxfinal += pow(2,i)*bcm_CSRMatrixNumNonzeros(A_array[i]); 
          }  
          wcmplxfinal=wcmplxfinal/bcm_CSRMatrixNumNonzeros(A);
          bcm_AMGHierarchyOpCmplxW(amg_hierarchy)=wcmplxfinal;
    } 

  time2=time_getWallclockSeconds()-time1;
  printf("Time for smoothed aggregation:  %e\n", time2);

  time1=time_getWallclockSeconds();
  if(coarse_solver == 9)
    {
      /* INSERT SLU HERE FOR THE MATRIX A_array[lev-1] */
#ifdef HAVE_SUPERLU
      factors_t   *LUfactors=bcm_AMGHierarchyLUfactors(amg_hierarchy);
      SuperMatrix *L, *U, B;
      int *perm_r; /* row permutations from partial pivoting */
      int *perm_c; /* column permutation vector */
      int *etree;  /* column elimination tree */
      SCformat *Lstore;
      NCformat *Ustore;
      double *b;
      int      i, panel_size, permc_spec, relax,nrhs=1, info;
      trans_t   trans = NOTRANS;
      mem_usage_t   mem_usage;
      superlu_options_t options;
      SuperLUStat_t stat;	      
      int  n       = bcm_CSRMatrixNumRows(A_array[lev-1]);
      int  nnz     = bcm_CSRMatrixNumNonzeros(A_array[lev-1]);
      if (LUfactors==NULL)
	{
	  SuperMatrix SA, AC;
	  GlobalLU_t Glu;   /* Not needed on return. */
	  double         *A_data  = bcm_CSRMatrixData(A_array[lev-1]);
	  int            *A_i     = bcm_CSRMatrixI(A_array[lev-1]);
	  int            *A_j     = bcm_CSRMatrixJ(A_array[lev-1]);
	  double         *values;
	  int            *rowind,*colptr;
	  
	  
	  
	  
	  /* Set the default input options. */
	  set_default_options(&options);
	  /* Initialize the statistics variables. */
	  StatInit(&stat);
	  dCompRow_to_CompCol(n,  n, nnz,
			      A_data, A_j, A_i,
			      &values,&rowind,&colptr);
	  dCreate_CompCol_Matrix(&SA, n, n, nnz, values, rowind, colptr,
				 SLU_NC, SLU_D, SLU_GE);
	  L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
	  U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
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
		 NULL, 0, perm_c, perm_r, L, U, &Glu, &stat, &info);
#elif defined(SLU_VERSION_4)
	  dgstrf(&options, &AC, relax, panel_size, etree,
		 NULL, 0, perm_c, perm_r, L, U, &stat, &info);
#else
	  choke_on_me;
#endif
	  
	  if ( info == 0 ) {
	    Lstore = (SCformat *) L->Store;
	    Ustore = (NCformat *) U->Store;
	    dQuerySpace(L, U, &mem_usage);
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
	      dQuerySpace(L, U, &mem_usage);
	      printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
		     mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	    }
	  }
	  
	  /* Save the LU factors in the factors handle */
	  LUfactors = (factors_t*) SUPERLU_MALLOC(sizeof(factors_t));
	  LUfactors->L = L;
	  LUfactors->U = U;
	  LUfactors->perm_c = perm_c;
	  LUfactors->perm_r = perm_r;
	  bcm_AMGHierarchyLUfactors(amg_hierarchy)= LUfactors;
	  /* Free un-wanted storage */
	  SUPERLU_FREE(etree);
	  Destroy_SuperMatrix_Store(&SA);
	  Destroy_CompCol_Permuted(&AC);
	  free(values);
	  free(rowind);
	  free(colptr);
	  StatFree(&stat);
	}
#endif      
    }
  time2=time_getWallclockSeconds()-time1;
  printf("Time for LU factorization of the coarsest %e\n", time2);

  printf("Current Coarsening Info\n");
  printf("Number of levels:  %d\n", bcm_AMGHierarchyNumLevels(amg_hierarchy));
  for(i=0; i<lev; i++)
    {
      printf("Size at level %d:  %d ; nnz %d\n", i, bcm_CSRMatrixNumRows(A_array[i]),
	     bcm_CSRMatrixNumNonzeros(A_array[i]));
    }
  printf("Current cmplx for V-cycle %e \n", bcm_AMGHierarchyOpCmplx(amg_hierarchy));
  printf("Current cmplx for W-cycle %e \n", bcm_AMGHierarchyOpCmplxW(amg_hierarchy));
  printf("Average Coarsening Ratio %e \n", bcm_AMGHierarchyAvgCratio(amg_hierarchy));

  return amg_hierarchy;
}
