#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 1038

int main (void)
{
  int i, warn;

  double lambda[N], flux[N], corr_flux[N], err_flux[N];
  
  FILE *data, *ready;
  data = fopen("Suma.txt", "r");
  ready = fopen("Omega_Cen.dat" , "w");
  
  for(i=0; i<N; i++)
    {  
      fscanf(data, "%lf\t %lf\t\n", &lambda[i], &flux[i]);
      
	  corr_flux[i] = 1.216746*flux[i]*pow(10,17);
	  err_flux[i] = corr_flux[i]*0.1;
	  fprintf(ready,"%lf\t %lf\t %lf\t %d\n", lambda[i], corr_flux[i], err_flux[i], 0);    
    }
    
  fclose(data);
  fclose(ready); 
}