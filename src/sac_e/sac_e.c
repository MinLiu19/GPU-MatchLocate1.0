/********************************************************
*	Calculate SAC data total energy (amplitude^2) withing specified time window
*	and normalize to single data point
*	Usage:
*		sac_e t1 t2 sac_files ...
*	Modified from lsac2, by Zhigang Peng, Fri Aug  9 11:13:54 PDT 2002
********************************************************/

#include <stdio.h>
#include <math.h>
#include "sac.h"

int main(int argc, char **argv) {
  SACHEAD	hd;
  int		i,j,n1,n2;
  float		*ar;
  float		t1, t2, t, energy;

  if (argc < 4) {
     fprintf(stderr, "Usage: %s t1 t2 sac_files ...\n\
  Calculate SAC data total energy (amplitude^2) withing specified time window\n\
  and normalize to single data point.\n",argv[0]);
     return -1;
  }

  energy = 0;

  sscanf(argv[1],"%f",&t1);
  sscanf(argv[2],"%f",&t2);
  for (i=3;i<argc;i++) {
     if ((ar = read_sac(argv[i],&hd)) == NULL) continue;
     n1= (int) ((t1-hd.b)/hd.delta);if(n1<1) n1=1;
     n2= (int) ((t2-hd.b)/hd.delta);if(n2>hd.npts-2) n2=hd.npts-2;
     if (n1>n2) {
        fprintf(stderr,"no time window for %s\n",argv[i]);
	continue;
     }
     t = n2-n1+1.;
     for(j=n1;j<n2;j++) {
	energy +=ar[j]*ar[j];	
     }
     energy = energy/t;
     printf("%s %e\n", argv[i],energy);
  }

  return 0;
}
