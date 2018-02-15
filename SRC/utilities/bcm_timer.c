/* 
                BootCMatch 
     Bootstrap AMG based on Compatible weighted Matching, version 0.9
    (C) Copyright 2017
                       Pasqua D'Ambra         IAC-CNR, IT
                       Panayot S. Vassilevki  Portland State University, OR USA
 
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
#include <time.h>
#include <unistd.h>
#include <sys/times.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>


double time_getWallclockSeconds(void)
{
#if 0
   struct tms usage;
   int wallclock = times(&usage);
   return(((double) wallclock)/((double) sysconf(_SC_CLK_TCK)));
#else
  struct timeval tt;
  struct timezone tz;
  double temp;
  if (gettimeofday(&tt,&tz) != 0) {
    fprintf(stderr,"Fatal error for gettimeofday ??? \n");
    exit(-1);
  }
  temp = ((double)tt.tv_sec) + ((double)tt.tv_usec)*1.0e-6;
  return(temp);
#endif   
}

double time_getCPUSeconds(void)
{
   clock_t cpuclock = clock();
   return(((double) (cpuclock))/((double) CLOCKS_PER_SEC));
}

double time_get_wallclock_seconds_(void)
{
   return(time_getWallclockSeconds());
}

double time_get_cpu_seconds_(void)
{
   return(time_getCPUSeconds());
}
