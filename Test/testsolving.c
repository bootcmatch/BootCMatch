/* 
                BootCMatch 
     Bootstrap AMG based on Compatible weighted Matching version 0.9
    (C) Copyright 2017
                       Pasqua D'Ambra    IAC-CNR, IT
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

#include "bcm.h"
#include "ioutil.h"

typedef struct {
  int matrixformat;
  char *matrixfile;
  char *rhsfile;
  char *solfile;
  int solver_type;
  int max_hrc;
  double conv_ratio;
  int matchtype;
  int aggrsweeps;
  int aggrtype;
  int max_levels;
  int cycle_type;
  int coarse_solver;
  int relax_type;
  int prerelax_sweeps;
  int postrelax_sweeps;
  int itnlim;
  double rtol;
} parms_t;

int get_inp_data(FILE *fp, parms_t *inparms)
{
  if (fp != NULL) {
    inparms->matrixformat      = int_get_fp(fp);
    inparms->matrixfile        = string_get_fp(fp);
    inparms->rhsfile           = string_get_fp(fp);
    if (strcmp(inparms->rhsfile,"NONE")==0) {
      free(inparms->rhsfile);
      inparms->rhsfile=NULL;
    }
    inparms->solfile           = string_get_fp(fp);
    if (strcmp(inparms->solfile,"NONE")==0) {
      free(inparms->solfile);
      inparms->solfile=NULL;
    }
    inparms->solver_type       = int_get_fp(fp);
    inparms->max_hrc           = int_get_fp(fp);
    inparms->conv_ratio        = double_get_fp(fp);
    inparms->matchtype         = int_get_fp(fp);
    inparms->aggrsweeps        = int_get_fp(fp);
    inparms->aggrtype          = int_get_fp(fp);
    inparms->max_levels        = int_get_fp(fp);
    inparms->cycle_type        = int_get_fp(fp);
    inparms->coarse_solver     = int_get_fp(fp);
    inparms->relax_type        = int_get_fp(fp);
    inparms->prerelax_sweeps   = int_get_fp(fp);
    inparms->postrelax_sweeps  = int_get_fp(fp);    
    inparms->itnlim            = int_get_fp(fp);
    inparms->rtol              = double_get_fp(fp);

  } else {
    return(-1);
  }
  return(0);
}

/* This program builds and applies a
   Bootstrap AMG with a desired convergence rate */

int  main(int argc, char *argv[])
{
   bcm_BootAMGBuildData *bootamg_data;
   bcm_AMGApplyData *amg_cycle;
   bcm_CSRMatrix *A;
   bcm_Vector *w;
   bcm_Vector *rhs, *Sol;
   bcm_Vector *soltrue;
   int i, *num_grid_sweeps, j, k;
   int size_rhs;
   bcm_Vector *res;
   double ratio, normnew, normold;
   parms_t inparms;
   FILE    *fp=NULL;

   fprintf(stdout,
	   "Welcome to BootCMatch version %s\n This is the testbootsolver program\n",
	   BCM_VERSION_STRING);

   if (argc > 1) {
     fp = fopen(argv[1], "r");
   } else {
     fp = stdin;
   }

   if (get_inp_data(fp,&inparms) != 0) {
     fprintf(stderr,"Error getting input parms \n");
     exit(1);
   }
   /* read sparse matrix from file */
   switch(inparms.matrixformat)
   {
     case 0:
   {
     fprintf(stderr,"Reading in COO\n");
     A=bcm_COO2CSRMatrixRead(inparms.matrixfile);
   }
   break;
     case 1:
   {
     fprintf(stderr,"Reading in MM\n");
     A=bcm_MM2CSRMatrixRead(inparms.matrixfile);
   }
   break;
     case 2:
   {
     fprintf(stderr,"Reading in CSR\n");
     A=bcm_CSRMatrixRead(inparms.matrixfile);
   }
   break;
   }
   int num_rows = bcm_CSRMatrixNumRows( A );
   int num_cols = bcm_CSRMatrixNumCols( A );
   int nnz = bcm_CSRMatrixNumNonzeros( A );
   printf("numrows: %d\n",num_rows);
   printf("numcols: %d\n",num_cols);
   printf("nnz: %d\n",nnz);

   /* read rhs vector from file */
   if (inparms.rhsfile != NULL) {
     if (inparms.matrixformat==1) rhs=bcm_VectorMMRead(inparms.rhsfile); 
     else rhs=bcm_VectorRead(inparms.rhsfile);
     size_rhs= bcm_VectorSize(rhs);
     printf("size_rhs: %d\n",size_rhs);
   } else {
     /* generate rhs vector */
     rhs=bcm_VectorCreate(num_rows);
     bcm_VectorInitialize(rhs);
     size_rhs=bcm_VectorSize(rhs);
     printf("size_rhs: %d\n",size_rhs);
     bcm_VectorSetConstantValues(rhs,1.0);
     //bcm_VectorSetRandomValues(rhs,1.0);
   }

 /* read reference solution vector from file */
   if (inparms.solfile != NULL) {
     if (inparms.matrixformat==1) soltrue=bcm_VectorMMRead(inparms.solfile);
     else soltrue=bcm_VectorRead(inparms.solfile);
     int size_soltrue= bcm_VectorSize(soltrue);
     printf("size_soltrue: %d\n",size_soltrue);
   }

   assert(num_rows == size_rhs); 

   /* initialize data structure for AMG building. See the setup
      routines for changing default values */
   bootamg_data = bcm_BootAMGBuildDataInitialize();
   bcm_AMGBuildData *amg_data;
   amg_data=bcm_BootAMGBuildDataCompData(bootamg_data);
   amg_cycle= bcm_AMGCycleInitialize();

   bcm_AMGBuildDataCSRMatrix(amg_data)=A;

   /* reading AMG algorithmic parameters */

   /* composition type  and maximum number of
      components for the bootstrap AMG */
   /* Single Hierarchy Type: matchingtype, number of aggregation sweeps
      for aggressive coarsening and aggregation type (smoothed or unsmoothed) */

   bcm_BootAMGBuildSetSolverType(bootamg_data, inparms.solver_type);
   bcm_BootAMGBuildSetMaxHrc(bootamg_data, inparms.max_hrc);
   bcm_BootAMGBuildSetDesiredRatio(bootamg_data, inparms.conv_ratio);
   bcm_AMGBuildSetAggMatchType(amg_data, inparms.matchtype);
   bcm_AMGBuildSetSweepNumber(amg_data, inparms.aggrsweeps);
   bcm_AMGBuildSetAggInterpType(amg_data, inparms.aggrtype);
   bcm_AMGBuildSetMaxLevels(amg_data, inparms.max_levels);
   bcm_AMGBuildSetCoarseSolver(amg_data, inparms.coarse_solver);
   bcm_AMGSetCycleType(amg_cycle, inparms.cycle_type);
   bcm_AMGSetRelaxType(amg_cycle, inparms.relax_type);
   bcm_AMGSetPreRelaxSteps(amg_cycle, inparms.prerelax_sweeps);
   bcm_AMGSetPostRelaxSteps(amg_cycle, inparms.postrelax_sweeps);


   printf("Composite Solver Type: %d\n",bcm_BootAMGBuildDataCompType(bootamg_data));
   printf("Max number of components: %d\n",bcm_BootAMGBuildDataMaxHrc(bootamg_data));
   printf("Desired Conv. Ratio: %e\n",bcm_BootAMGBuildDataDesRatio(bootamg_data));
   printf("Matching Type: %d\n",bcm_AMGBuildDataAggMatchType(amg_data));
   printf("Aggregation sweeps: %d\n",bcm_AMGBuildDataSweepNumber(amg_data));
   printf("Aggregation Type: %d\n",bcm_AMGBuildDataAggInterpType(amg_data));
   printf("max_levels: %d\n",bcm_AMGBuildDataMaxLevels(amg_data));
   printf("coarse_solver: %d\n",bcm_AMGBuildDataCoarseSolver(amg_data));
   printf("cycle_type: %d\n",bcm_AMGApplyDataCycleType(amg_cycle));
   printf("relax_type: %d\n",bcm_AMGApplyDataRelaxType(amg_cycle));
   printf("prerelax_sweeps: %d\n",bcm_AMGApplyDataPreRelax(amg_cycle));
   printf("postrelax_sweeps: %d\n",bcm_AMGApplyDataPostRelax(amg_cycle));

   printf("itnlim: %d\n",inparms.itnlim);
   printf("rtol: %e\n",inparms.rtol);

   if (argc > 1) fclose(fp);

   /* set maximum coarse size */
   
   int maxcoarse=40*pow((double)num_rows,(double)1/3);
   bcm_AMGBuildSetMaxCoarseSize(amg_data, maxcoarse);
   printf("maxcoarsesize: %d\n",bcm_AMGBuildDataMaxCoarseSize(amg_data));

   /* No relaxation of smooth vectors in AMG building */
   int CRit=0;
   bcm_AMGBuildSetCRIterations(amg_data, CRit);
   printf("CRit: %d\n",bcm_AMGBuildDataCRIterations(amg_data));

   /* initialize num_grid_sweeps parameter w.r.t. the number
      of levels and the chosen cycle.
      In the final release we have to manage this in a setup routine after hierarchy building */
   
   num_grid_sweeps = (int *) calloc(inparms.max_levels-1, sizeof(int));
   for(i=0; i<inparms.max_levels-1; i++) num_grid_sweeps[i]=1;
   switch(inparms.cycle_type)  {
     case 1: /* H-cycle */
       {
	 for(i=0; i<inparms.max_levels-1; i++)  {
	   j=i%2; /*step is fixed to 2; it can be also different */
	   if(j==0) num_grid_sweeps[i]=2; /* if num_grid_sweeps is 2, we are using a hybrid V-W cycle */
	 }
       }
       break;
   case 2: /* W-cycle */
     {
       for(i=0; i<inparms.max_levels-2; i++) num_grid_sweeps[i]=2;
     }
     break;
   }
   bcm_AMGApplyDataGridSweeps(amg_cycle)=num_grid_sweeps;


   /* set arbitrary initial (smooth) vector: generally
      we use unitary vector or random vector */

   w=bcm_VectorCreate(num_rows);
   bcm_VectorInitialize(w);
   int w_size=bcm_VectorSize(w);
   printf("wsize: %d\n",w_size);
   /*bcm_VectorSetRandomValues(w,1.0);*/
   bcm_VectorSetConstantValues(w,1.0); 

   bcm_AMGBuildDataSmoothVector(amg_data)= w;

   /* start bootstrap process */

   printf("Bootstrap starting \n");

   bcm_BootAMG *boot_amg;
   double time1=time_getWallclockSeconds();
   boot_amg=bcm_Bootstrap(bootamg_data,amg_cycle);
   double time2=time_getWallclockSeconds()-time1;
   
   printf("Bootstrap ended\n");
   
   bcm_AMGHierarchy **Harray;
   
   printf("Number of components:  %d\n", bcm_BootAMGNHrc(boot_amg));
   printf("Estimated convergence %e \n", bcm_BootAMGEstRatio(boot_amg));
   printf("Information on the Components\n");
   Harray=bcm_BootAMGHarray(boot_amg);
   for(k=0;k<bcm_BootAMGNHrc(boot_amg); k++)  {
     printf("Component:  %d\n", k);
     printf("Number of levels:  %d\n", bcm_AMGHierarchyNumLevels(Harray[k]));
     printf("Operator cmplx for V-cycle %e \n", bcm_AMGHierarchyOpCmplx(Harray[k]));
     printf("Operator cmplx for W-cycle %e \n", bcm_AMGHierarchyOpCmplxW(Harray[k]));
   }
   printf("Wall Clock Time for Building:  %e\n", time2);
   /* printing statistics */
   double avgcratio=0.0;
   double avgcmpx=0.0;
   double avgwcmpx=0.0;
   double avgnumlev=0.0;
   for (k=0;k<bcm_BootAMGNHrc(boot_amg); k++) {
     avgnumlev=avgnumlev+bcm_AMGHierarchyNumLevels(Harray[k]);
     avgcmpx=avgcmpx+bcm_AMGHierarchyOpCmplx(Harray[k]);
     avgwcmpx=avgwcmpx+bcm_AMGHierarchyOpCmplxW(Harray[k]);
     avgcratio=avgcratio+bcm_AMGHierarchyAvgCratio(Harray[k]);
   }
   avgnumlev=avgnumlev/bcm_BootAMGNHrc(boot_amg);
   avgcmpx=avgcmpx/bcm_BootAMGNHrc(boot_amg);
   avgwcmpx=avgwcmpx/bcm_BootAMGNHrc(boot_amg);
   avgcratio=avgcratio/bcm_BootAMGNHrc(boot_amg);
   printf("Average of number of levels:  %e\n", avgnumlev);
   printf("Average of operator complexity for V-cycle:  %e\n", avgcmpx);
   printf("Average of operator complexity for W-cycle:  %e\n", avgwcmpx);
   printf("Average of coarsening ratio:  %e\n", avgcratio);
   
   /* generate initial vector */
   Sol=bcm_VectorCreate(num_rows);
   bcm_VectorInitialize(Sol);
   bcm_VectorSetConstantValues(Sol,0.0);
   
   printf("starting composite amg solver\n");
   int ierr;
   time1=time_getWallclockSeconds();
   ierr=bcm_boot_solver(bootamg_data, boot_amg, amg_cycle,rhs,
			   Sol,inparms.itnlim,inparms.rtol);
   time2=time_getWallclockSeconds()-time1;
   printf("Wall Clock Time for compamg_solver:  %e\n", time2);
   
   char *file_out="sol.mtx";
   bcm_VectorPrint(Sol,file_out);

     if (inparms.solfile != NULL)        {
     double alpha=-1.0;
     bcm_Vector *diff;
     diff=bcm_VectorCreate(num_rows);
     bcm_VectorInitialize(diff);

     bcm_VectorCopy(Sol,diff);
     bcm_VectorAxpy(alpha,soltrue,diff);
     double norm2diff;
     norm2diff=bcm_VectorNorm(diff);
     double norm2reldiff;
     norm2reldiff=norm2diff/bcm_VectorNorm(soltrue);
     printf("Solution error (Euclidean norm):  %e\n", norm2diff);
     printf("Solution relative error (Euclidean norm):  %e\n", norm2reldiff);

     bcm_VectorDestroy(diff);
     bcm_VectorDestroy(soltrue);

     }

   
   bcm_BootAMGBuildDataDestroy(bootamg_data);
   bcm_AMGApplyDataDestroy(amg_cycle);
   bcm_BootAMGDestroy(boot_amg);
   bcm_VectorDestroy(rhs);
   bcm_VectorDestroy(Sol);
   bcm_VectorDestroy(w);
   bcm_CSRMatrixDestroy(A);
   free(num_grid_sweeps);
   return (0);
}
