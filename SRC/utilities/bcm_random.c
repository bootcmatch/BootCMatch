/* 
                BootCMatch 
     Bootstrap AMG based on Compatible Matching version 0.9
    (C) Copyright 2017
                       Pasqua D'Ambra         IAC-CNR, IT
                       Salvatore Filippone    Cranfield University, UK
                       Panayot S. Vassilevki CASC-LLNL
 
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
 * This file contains routines that implement a pseudo-random number generator
 * detailed in the following paper.
 *
 * @article{RNG_Park_Miller,
 *   author = {S. K. Park and K. W. Miller},
 *   title = {Random number generators: good ones are hard to find},
 *   journal = {Commun. ACM},
 *   volume = {31},
 *   number = {10},
 *   year = {1988},
 *   pages = {1192--1201},
 * }
 *
 * This RNG has been shown to appear fairly random, it is a full period
 * generating function (the sequence uses all of the values available to it up
 * to 2147483647), and can be implemented on any architecture using 32-bit
 * integers. The implementation in this file will not overflow for 32-bit
 * arithmetic, which all modern computers should support.
 *
 * @author David Alber
 * @date March 2005
 *
 *****************************************************************************/

/*--------------------------------------------------------------------------
 * Static variables
 *--------------------------------------------------------------------------*/

static int Seed = 13579;

#define a  16807
#define m  2147483647
#define q  127773
#define r  2836

/*--------------------------------------------------------------------------
 * Initializes the pseudo-random number generator to a place in the sequence.
 *
 * @param seed an int containing the seed for the RNG.
 *--------------------------------------------------------------------------*/

void  bcm_SeedRand( int seed )
{
   Seed = seed;
}

/*--------------------------------------------------------------------------
 * Computes the next pseudo-random number in the sequence using the global
 * variable Seed.
 *
 * @return a double containing the next number in the sequence divided by
 * 2147483647 so that the numbers are in (0, 1].
 *--------------------------------------------------------------------------*/

double  bcm_Rand()
{
   int  low, high, test;

   high = Seed / q;
   low = Seed % q;
   test = a * low - r * high;
   if(test > 0)
   {
      Seed = test;
   }
   else
   {
      Seed = test + m;
   }

   return ((double)(Seed) / m);
}
