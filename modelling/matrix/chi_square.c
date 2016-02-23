#include <stdio.h>
#include <math.h>

#define N 105

int main(){
  
  FILE *mod,*sigma;
  int warn, i, j;
  double R[12], sig_p[12], Re[N], sig_pe[N], shi2, shi_t;
  
  sigma=fopen("sigma.dat","r"); 
  mod=fopen("model.dat","r");
  
  for(j=0.0; j<12; j++)
    {        
      fscanf(sigma, "%lf %lf", &R[j], &sig_p[j]);
    }      
  
  fclose(sigma);
  
    for(i=0.0; i<N; i++)
    {        
      fscanf(mod, "%lf %lf", &Re[i], &sig_pe[i]);
      
        if( Re[i] < 0.1 && sig_pe[i] < 0.1 )
	    {
	      printf("%lf\t %lf\n", Re[i], sig_pe[i]);
	       //printf("%lf\n", shi_t);       
	    }
    }      
  
  fclose(mod);
  
  for(j=0.0; j<12.0; j=j+1.0)
    {
      for(i = 0.0;i<N; i++)
	{
	  	 
	    if( R[j] - Re[i] < 0.2 && R[j] - Re[i] > -0.2 )
	    {   
	      //printf("%lf\t %lf\t %lf\t %lf\n", R[j], sig_p[j], Re[i], sig_pe[i]);
	       shi2 = (sig_p[j]-sig_pe[i])*(sig_p[j]-sig_pe[i]);
	       shi_t = shi_t + shi2;
	       //printf("%lf\n", shi2);
	       //printf("%lf\n", shi_t);       
	    }
	}
    }
    printf("%lf\n", shi_t);
 
  return(warn);
}
