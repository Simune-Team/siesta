/* $Id: glib.c,v 1.2 1999/01/31 12:05:51 emilio Exp $ */
/*
  Auxiliary routines
*/

/*
   From elgin@claudia.spectral.com (Jim Elgin)
*/

/*  supply "etime" fortran routine

 NAME
      etime - return elapsed execution time

 SYNOPSIS
      REAL function etime (tarray)
      REAL tarray(2)

 DESCRIPTION
      This routine returns elapsed runtime in seconds for the calling
      process.

      The argument array returns user time in the first element and system
      time in the second element.  The function value is the sum of user and
      system time.

      The resolution of all timing is 1/CLK_TCK seconds, where CLK_TCK is
      processor dependent.
*/
#include <unistd.h>
#include <sys/times.h>

float etime_(tarray)
float *tarray;
{
struct tms buf;
float t1, t2, den, tot;

times(&buf);
t1 = buf.tms_utime;
t2 = buf.tms_stime;
den = sysconf(_SC_CLK_TCK);
*tarray = t1/den;
*(tarray+1) = t2/den;
tot = *tarray + *(tarray+1);
return tot;
}
#include<time.h>
void fdate_(utime)
char utime[24];
{
int i;
time_t t;
t=time (NULL);
for (i=0; i<24; i++)
 utime[i]= *(ctime(&t)+i);
}

#include <stdlib.h>
void exit_(status)
int *status;
{
  exit(*status);
}

