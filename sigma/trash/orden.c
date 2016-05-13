#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#define N 3749
#define N_bin 250

int main (void)
{
  int i, j, warn, counter;
  double lim[11];
  counter = 0.0;
  
  double pos[N],vel[N],err[N],vel_sq, sum;
  double mean, median, skew, sd, kurtosis, error, sum1, sum2, temp1=0, temp2=0, temp3=0;
  
  FILE *velo, *sig, *script, *order;
  velo = fopen("velocity.dat", "r");
  sig = fopen("sigma.dat","w");
  order = fopen("organized.dat","w");
  
  for(i=0; i<N; i++)
    {        
      fscanf(velo, "%lf %lf %lf", &pos[i], &vel[i], &err[i]);
    }
  
  for(i=0; i<(N-1); ++i)
    {
      for(j=0; j<N-1-i; ++j)
	{
	  if(pos[j] < pos[j+1])
	    {		
	      temp1 = pos[j+1];  
	      pos[j+1]=pos[j];          
	      pos[j]=temp1; 
	      
	      temp2=vel[j+1];
	      vel[j+1]=vel[j];
	      vel[j]=temp2;
	      
	      temp3=err[j+1];
	      err[j+1]=err[j];
	      err[j]=temp3;	    
	    }	  	  
	}
      fprintf(order, "%lf %lf %lf\n", pos[j], vel[j], err[j]);
    }
  fclose(velo);      
  fclose(order);
  return(warn);
}
