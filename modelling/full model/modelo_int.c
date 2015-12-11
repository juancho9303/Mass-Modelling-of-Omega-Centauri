#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#define N 100
#define pi 3.141592
#define M 1.989
#define gamma 1.2
#define G 43007.1

double X_s,I_R,s,R,a,sig_p,beta,rho,re1,re2;

double int_1 (double r, void * params) {
  //  double a = *(double *) params;
  double int_1 = 1.0/(r*pow((r+a),5));
  return int_1;
}

void main()
{
  
  FILE *mod,*script;
  
  mod=fopen("model.dat","w");
    
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (2000);
    
      double result1, error1;

  gsl_function F1;
  F1.function = &int_1;
  //  F1.params = &a;
  
  gsl_integration_qags (&F1, 1, 10000, 0, 1e-4, 2000,
                        w, &result1, &error1); 
  
  printf ("result          = % .18f\n", result1);
  printf ("estimated error = % .18f\n", error1);
  printf ("intervals       = %zu\n", w->size);
  
  gsl_integration_workspace_free (w);
  
  for (R = 0.01; R < N; R=R+0.01)
    {      
      re1 = result1;
      
      double int_2 (double r, void * params) {
	//  double a = *(double *) params;
	double int_2 = (1.0 - beta*((R*R)/(r*r)))*(re1)*((r)/(sqrt(r*r-R*R)));
	//  return int_2;
      } 
      
      if (R < 1.0)
	{
	  
	  gsl_integration_workspace * z 
	    = gsl_integration_workspace_alloc (2000);
	    
	  double result2, error2;
	  
	  gsl_function F2;
	  F2.function = &int_2;
	  //  F2.params = &a;
	  
	  gsl_integration_qags (&F2, R, 100000, 0, 1e-4, 2000,
				z, &result2, &error2); 
	  
	  printf ("result = % .18f\n", result2);
	  
	  gsl_integration_workspace_free (z);
	  
	  re2 = result2;
	  beta = 2.0;
	  a = 2.0;
	  s = R/a;
	  X_s = (1.0 / sqrt(1.0 - s*s)) * log((1.0+(sqrt(1-s*s)))/(s));
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
	  
	  double result2, error2;
	  
	  gsl_function F2;
	  F2.function = &int_2;
	  //  F2.params = &a;
	  
	  gsl_integration_qags (&F2, R, 100000, 0, 1e-4, 2000,
				z, &result2, &error2); 
	  
	  printf ("result  = % .18f\n", result2);
	  
	  gsl_integration_workspace_free (z);
	  
	  re2 = result2;
	  beta = 2.0;
	  a = 2.0;
	  s = R/a;
	  X_s = (1.0 / sqrt(s*s - 1.0)) * acos(1/s);
	  I_R = (M/(2.0*pi*a*a*gamma*(1.0-s*s)*(1.0-s*s)))*((2.0+s*s)*X_s-3.0);
	  sig_p = ((2.0*G*M*M*a)/(gamma*I_R*2.0*pi))*re2;
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
