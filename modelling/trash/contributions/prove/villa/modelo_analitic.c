#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#define N 10
#define M 2
#define GAMMA 1.2
#define G 43007.1

struct param{
  double a;
  double a_ste;
  double beta;
  double M_ste;
  double a_dm;
  double M_dm;
};

double R;

double A_r (double r, void * params) 
{
  struct param parameters = *(struct param *) params;
  double a  = parameters.a;
  double beta = parameters.beta;  
  double a_ste = parameters.a_ste;
  double M_ste = parameters.M_ste;
  double A_r = (M_ste*M_ste*a_ste)/(r*pow((r+a_ste),5));
  return A_r;
}

double B_r (double r, void * params) 
{
  struct param parameters = *(struct param *) params;
  double a  = parameters.a;
  double beta = parameters.beta;  
  double a_ste = parameters.a_ste;
  double M_ste = parameters.M_ste;
  double a_dm = parameters.a_dm;
  double M_dm = parameters.M_dm;
  double B_r = (M_ste*M_dm*a_ste)/(r*pow((r+a_ste),3)*pow((r+a_dm),2));
  return B_r;
}

double C_r (double r, void * params) 
{
  struct param parameters = *(struct param *) params;
  double a  = parameters.a;
  double beta = parameters.beta;  
  double a_ste = parameters.a_ste;
  double M_ste = parameters.M_ste;
  double a_dm = parameters.a_dm;
  double M_dm = parameters.M_dm;
  double C_r = (M_dm*M_ste*a_dm)/(r*pow((r+a_dm),3)*pow((r+a_ste),2));
  return C_r;
}

double D_r (double r, void * params) 
{
  struct param parameters = *(struct param *) params;
  double a  = parameters.a;
  double beta = parameters.beta;  
  double a_dm = parameters.a_dm;
  double M_dm = parameters.M_dm;
  double D_r = (M_dm*M_dm*a_dm)/(r*pow((r+a_dm),5));
  return D_r;
}

int main(){
  
  FILE *mod,*script;
  int warn;
  double X_s,I_R,s,R,sig_p,rho,re1,re2,re3,re4,re5,result1,error1,result2,error2,result3,error3,result4,error4,result5,error5;
  int Nint = 2000; // Numero de intervalos
  struct param params; // structura de parametros
  gsl_integration_workspace *z = gsl_integration_workspace_alloc(Nint);
  gsl_function F1;
  gsl_function F2;
  gsl_function F3;
  gsl_function F4;
  gsl_function F5;
  
  // parametros del problema (o algo asi)
  double a = 1.0;
  double a_ste = 1.0;
  double beta = 2.0;F1.function = &A_r;
    F1.params = &params;

    gsl_integration_qags(&F1, 1, 1000, 0, 1e-4, Nint, z, &result1, &error1); 
    //printf ("result  = % .18f\n", result1);
  double M_ste = 2.0;
  double a_dm = 1.0;
  double M_dm = 2.0;

  mod=fopen("model.dat","w");
  params.a = a;
  params.a_ste = a_ste;
  params.beta = beta;
  params.M_ste = M_ste;
  params.a_dm = a_dm;
  params.M_dm = M_dm;
      
  for(R = 0.01; R < N; R=R+0.01){

    F1.function = &A_r;
    F1.params = &params;

    gsl_integration_qags(&F1, 1, 1000, 0, 1e-4, Nint, z, &result1, &error1); 
    //printf ("result  = % .18f\n", result1);
    
    F2.function = &B_r;
    F2.params = &params;

    gsl_integration_qags(&F2, 1, 1000, 0, 1e-4, Nint, z, &result2, &error2); 
    //printf ("result  = % .18f\n", result2);
    
    F3.function = &C_r;
    F3.params = &params;

    gsl_integration_qags(&F3, 1, 1000, 0, 1e-4, Nint, z, &result3, &error3); 
    //printf ("result  = % .18f\n", result2);
    
    F4.function = &D_r;
    F4.params = &params;

    gsl_integration_qags(&F4, 1, 1000, 0, 1e-4, Nint, z, &result4, &error4); 
    //printf ("result  = % .18f\n", result2);

    double tot_r (double r, void * params) 
{
  struct param parameters = *(struct param *) params;
  double a  = parameters.a;
  double beta = parameters.beta;  
  double a_dm = parameters.a_dm;
  double M_dm = parameters.M_dm;  
  double tot_r = (1-beta*R*R/r*r)*(r/(sqrt(r*r-R*R)))*(result1+result2+result3+result4);
  return tot_r;
}
    
    F5.function = &tot_r;
    F5.params = &params;

    gsl_integration_qags(&F5, R, 1000, 0, 1e-4, Nint, z, &result5, &error5); 
    //printf ("result  = % .18f\n", result5);
    
    if (R < 1.0)
      {      
	
	re5 = result5;
	s = R/a;
	X_s = (1.0 / sqrt(1.0 - s*s)) * log((1.0+(sqrt(1.0-s*s)))/(s));
	I_R = (M/(2.0*M_PI*a*a*GAMMA*(1.0-s*s)*(1.0-s*s)))*((2.0+s*s)*X_s-3.0);
	sig_p = ((2.0*G*M*M*a)/(GAMMA*I_R*2.0*M_PI))*re5;
	//printf("sig_p^2(R) %6lf\n", sig_p);
	fprintf(mod,"%6lf\t %6lf\n", s, sig_p);
      }
    else if(R > 1.01)
      {	       	  
	re5 = result5;
	s = R/a;
	X_s = (1.0 / sqrt(s*s - 1.0)) * acos(1.0/s);
	I_R = (M/(2.0*M_PI*a*a*GAMMA*(1.0-s*s)*(1.0-s*s)))*((2.0+s*s)*X_s-3.0);
	sig_p = ((2.0*G*M*M*a)/(GAMMA*(I_R)*2.0*M_PI))*re5;
	//printf("sig_p^2(R) %6lf\n", sig_p);
	fprintf(mod,"%6lf\t %6lf\n", s, sig_p);
      }      
  }

  gsl_integration_workspace_free(z);
  
  fclose(mod);
  
  script = fopen( "script.gpl", "w" );
  fprintf(script, "set terminal png\nset output 'sigma.png'\nset nokey\n");
  fprintf( script, "set title 'Sigma Proyectada vs R'\n" );
  fprintf( script, "set xlabel 'Radio (Kpc)'\n" );
  fprintf( script, "set ylabel 'Radio (Sigma_p (R))'\n" );
  fprintf( script, "set grid\n" );
  fprintf( script, "plot 'model.dat' u 1:2 w d\n");
  fclose(script);
  warn = system("gnuplot script.gpl");
  
  return(warn);
}