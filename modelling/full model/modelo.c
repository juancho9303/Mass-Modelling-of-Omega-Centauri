#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#define N 10
#define pi 3.141592
#define M 5000
#define gamma 1.2
#define G 6.67

double X_s,I_R,s,R,a,sig_p,beta,rho;

void main()
{
  
//double X_s,I_R,s,R,a,sig_p,beta,rho;
  
  FILE *iso,*script;
  
  iso=fopen("isotropico.dat","w");
  
  for (R = 0.01; R < N; R=R+0.01)
    {
      
      if (R < 1.0)
	{
	  beta = 2.0;
	  a = 2.0;
	  s = R/a;
	  X_s = (1.0 / sqrt(1.0 - s*s)) * log((1.0+(sqrt(1-s*s)))/(s));
	  I_R = (M/(2.0*pi*a*a*gamma*(1.0-s*s)*(1.0-s*s)))*((2.0+s*s)*X_s-3.0);
	  sig_p = (2.0*G*M*M*a)/(gamma*I_R*2.0*pi);
	  printf("sig_p^2(R) %6lf\n", sig_p);
	  fprintf(iso,"%6lf\t %6lf\n", s, sig_p);
	}
      else if(R > 1.01)
	{
	  beta = 2.0;
	  a = 2.0;
	  s = R/a;
	  X_s = (1.0 / sqrt(s*s - 1.0)) * acos(1/s);
	  I_R = (M/(2.0*pi*a*a*gamma*(1.0-s*s)*(1.0-s*s)))*((2.0+s*s)*X_s-3.0);
	  sig_p = (2.0*G*M*M*a)/(gamma*I_R*2.0*pi);
	  printf("sig_p^2(R) %6lf\n", sig_p);
	  fprintf(iso,"%6lf\t %6lf\n", s, sig_p);
	}      
    }
  
  fclose(iso);
  
  script = fopen( "script.gpl", "w" );
  fprintf(script, "set terminal png\nset output 'sigma.png'\nset nokey\n");
  fprintf( script, "set title 'Sigma Proyectada vs R'\n" );
  fprintf( script, "plot 'isotropico.dat' u 1:2 w d\n");
  fclose(script);
  system("gnuplot script.gpl");
  
}