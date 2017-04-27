/* 
                BootCMatch 
     Bootstrap AMG based on Compatible Matching version 0.9
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
 * functions for BootAMG data structures
 *
 *****************************************************************************/

#include "bcm_bootamg.h"

/*--------------------------------------------------------------------------
 * functions for bcm_BootAMGBuildData
 *--------------------------------------------------------------------------*/

void *
bcm_BootAMGBuildDataInitialize()
{
   bcm_BootAMGBuildData *bootamg_data;

   /* setup params */
   int      max_hrc; /* max number of hierarchies to be built*/
   double   conv_ratio; /* desired convergence ratio */
   double   estimated_ratio; /* estimated convergence ratio after built */
   int      solver_type; /* type of composition for applying composite solver */
  /* solver_type =0 -> multiplicative composition
   * solver_type =1 -> symmetrized multiplicative comsposition
   * solver_type =2 -> additive composition */
   int      solver_it; /* number of iterations to be applied for conv. ratio estimating */

   /*-----------------------------------------------------------------------
    * Setup default values for parameters
    *-----------------------------------------------------------------------*/

   /* setup params */
     max_hrc=10;
     conv_ratio=0.80; 
     solver_type=1; 
     solver_it=15; 

   /*-----------------------------------------------------------------------
    * Create the bcm_BootAMGBuildData structure and return
    *-----------------------------------------------------------------------*/

   bootamg_data = (bcm_BootAMGBuildData *) calloc(1, sizeof(bcm_BootAMGBuildData));

   bcm_BootAMGBuildSetMaxHrc(bootamg_data, max_hrc);
   bcm_BootAMGBuildSetDesiredRatio(bootamg_data, conv_ratio);
   bcm_BootAMGBuildSetSolverType(bootamg_data, solver_type);
   bcm_BootAMGBuildSetSolverIt(bootamg_data, solver_it);

   bcm_AMGBuildData *amg_data;
   amg_data=bcm_AMGInitialize();
   bcm_BootAMGBuildDataCompData(bootamg_data)=amg_data;

   return (void *) bootamg_data;
}
/*--------------------------------------------------------------------------
 * Routines to set the setup parameters
 *--------------------------------------------------------------------------*/

int
bcm_BootAMGBuildSetMaxHrc( void *data,
                       int   max_hrc )
{
   int ierr = 0;
   bcm_BootAMGBuildData  *bootamg_data = data;
 
   bcm_BootAMGBuildDataMaxHrc(bootamg_data) = max_hrc;

   return (ierr);
}

int
bcm_BootAMGBuildSetDesiredRatio( void     *data,
                        double       conv_ratio )
{
   int ierr = 0;
   bcm_BootAMGBuildData  *bootamg_data = data;

   bcm_BootAMGBuildDataDesRatio(bootamg_data) = conv_ratio;

   return (ierr);
}

int
bcm_BootAMGBuildSetSolverType( void     *data,
                     int       solver_type )
{
   int ierr = 0;
   bcm_BootAMGBuildData  *bootamg_data = data;
 
   bcm_BootAMGBuildDataCompType(bootamg_data) = solver_type;

   return (ierr);
} 

int
bcm_BootAMGBuildSetSolverIt( void     *data,
                           int      solver_it )
{
   int ierr = 0;
   bcm_BootAMGBuildData  *bootamg_data = data;

   bcm_BootAMGBuildDataCompIt(bootamg_data) = solver_it;

   return (ierr);
}


int
bcm_BootAMGBuildDataDestroy(void *data)

{
   int ierr=0;
   bcm_BootAMGBuildData  *bootamg_data = data;

   if(bootamg_data)
   {
     bcm_AMGBuildDataDestroy(bcm_BootAMGBuildDataCompData(bootamg_data));
     free(bootamg_data);
   }
   return ierr;
}

/*--------------------------------------------------------------------------
 * functions for bcm_BootAMG
 *--------------------------------------------------------------------------*/
void *
bcm_BootAMGCreate(int max_hrc)

{
  bcm_BootAMG *boot_amg;

  boot_amg = (bcm_BootAMG *) calloc(1, sizeof(bcm_BootAMG));

  bcm_BootAMGHarray(boot_amg)=NULL;
  bcm_BootAMGNHrc(boot_amg)=max_hrc;

  return boot_amg;
}

int 
bcm_BootAMGHarrayDestroy(int nhrc, int start, bcm_AMGHierarchy **Harray)
{
   int ierr=0;
   if(Harray) {
     int ih;
     for (ih=start; ih<nhrc; ih++) {       
       bcm_AMGHierarchyDestroy(Harray[ih]);
     }
   }
   return ierr;
}

int 
bcm_BootAMGDestroy(bcm_BootAMG *boot_amg)

{
   int ierr=0;

   if(boot_amg)
   {
     bcm_BootAMGHarrayDestroy(bcm_BootAMGNHrc(boot_amg),0,bcm_BootAMGHarray(boot_amg));
     free(bcm_BootAMGHarray(boot_amg));
     free(boot_amg);
   }
   return ierr;
}

int 
bcm_BootAMGInitialize(bcm_BootAMG *boot_amg)
{
   int n_hrc=bcm_BootAMGNHrc(boot_amg);

   int ierr=0;

   if(n_hrc)
   {
   if(! bcm_BootAMGHarray(boot_amg))
      bcm_BootAMGHarray(boot_amg)= (bcm_AMGHierarchy **) calloc(n_hrc, sizeof(bcm_AMGHierarchy));
   }
   bcm_BootAMGEstRatio(boot_amg)=1.0;
   return ierr;
}
