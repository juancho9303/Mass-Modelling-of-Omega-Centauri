#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#define N 20
#define pi 3.141592
#define M 10e-3
#define gamma 1.2
#define G 43007.1

double X_s,I_R,s,R,a,sig_p,beta,rho,re1,re2,result1,error1,result2,error2;

void main()
{
  
  FILE *mod,*script;
  
  mod=fopen("model.dat","w");
      
  for (R = 0.01; R < N; R=R+0.01)
    {      
      
      double int_2 (double r, void * params) {
	//  double a = *(double *) params;
	double int_2 = (1.0 - beta*((R*R)/(r*r)))*(((12.0*pow((a+r),4.0)*log((a+r)/(r)))-a*(25.0*a*a*a+52.0*a*a*r+42.0*a*r*r+12.0*r*r*r))/(12.0*pow(a,5)*pow((a+r),4)))*((r)/(sqrt(r*r-R*R)));
	  //return int_2;
      } 
      
      if (R < 1.0)
	{
	  
	  gsl_integration_workspace * z 
	    = gsl_integration_workspace_alloc (2000);
	  
	  gsl_function F2;
	  F2.function = &int_2;
	  //  F2.params = &a;
	  
	  gsl_integration_qags (&F2, 1, 1000, 0, 1e-4, 2000,
				z, &result2, &error2); 
	  
	  printf ("result = % .18f\n", result2);
	  
	  gsl_integration_workspace_free (z);
	  
	  re2 = result2;
	  beta = 2.0;
	  a = 2.0;
	  s = R/a;
	  X_s = (1.0 / sqrt(1.0 - s*s)) * log((1.0+(sqrt(1.0-s*s)))/(s));
	  I_R = (M/(2.0*pi*a*a*gamma*(1.0-s*s)*(1.0-s*s)))*((2.0+s*s)*X_s-3.0);
	  sig_p = ((2.0*G*M*M*a)/(gamma*I_R*2.0*pi))*re2;
	  //printf("sig_p^2(R) %6lf\n", sig_p);
	  fprintf(mod,"%6lf\t %6lf\n", s, sig_p);
	}
      else if(R > 1.01)
	{	  
	  //  int_2 = 1.0/(R*pow((R+a),3));
	  
	  gsl_integration_workspace * z 
	    = gsl_integration_workspace_alloc (2000);
	  
	  gsl_function F2;
	  F2.function = &int_2;
	  //  F2.params = &a;
	  
	  gsl_integration_qags (&F2, 1, 1000, 0, 1e-4, 2000,
				z, &result2, &error2); 
	  
	  printf ("result = % .18f\n", result2);
	  
	  gsl_integration_workspace_free (z);
	  
	  re2 = result2;
	  beta = 2.0;
	  a = 2.0;
	  s = R/a;
	  X_s = (1.0 / sqrt(s*s - 1.0)) * acos(1.0/s);
	  I_R = (M/(2.0*pi*a*a*gamma*(1.0-s*s)*(1.0-s*s)))*((2.0+s*s)*X_s-3.0);
	  sig_p = ((2.0*G*M*M*a)/(gamma*(I_R)*2.0*pi))*re2;
	  //printf("sig_p^2(R) %6lf\n", sig_p);
	  fprintf(mod,"%6lf\t %6lf\n", s, sig_p);
	}      
    }
  
  fclose(mod);
  
  script = fopen( "script.gpl", "w" );
  fprintf(script, "set terminal png\nset output 'sigma.png'\nset nokey\n");
  fprintf( script, "set title 'Sigma Proyectada vs R'\n" );
  fprintf( script, "plot 'model.dat' u 1:2 w d\n");
  fclose(script);
  system("gnuplot script.gpl");
  
}