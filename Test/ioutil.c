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
#include "ioutil.h"
#define BUFSIZE 2048
#define DELIM "%"

char linebuffer[BUFSIZE+1];

int    int_get_fp(FILE *fp)
{
  int temp;
  char *token;
  fgets(linebuffer,BUFSIZE,fp);
  token = strtok(linebuffer,DELIM);
  sscanf(token,"%d",&temp);
  return(temp);  
}

double    double_get_fp(FILE *fp)
{
  double temp;
  char *token;
  fgets(linebuffer,BUFSIZE,fp);
  token = strtok(linebuffer,DELIM);
  sscanf(token,"%lf",&temp);
  return(temp);  
}

char* string_get_fp(FILE *fp)
{
  char *token1, *token2;
  fgets(linebuffer,BUFSIZE,fp);
  token1 = strtok(linebuffer,DELIM);
  token2 = strtok(token1," ");
  return(strdup(token2));  
}
