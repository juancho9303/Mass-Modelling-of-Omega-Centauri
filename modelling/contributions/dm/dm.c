#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#define N 20
//#define M 10e-3
#define GAMMA 1.2
#define G 43007.1

struct param{
  //double a;
  double a_dm;
  double M_dm;
  double a_ste;
  double M_ste;
  double beta;
};

double R;

//Defino la funcion a integrar de la contribucion de materia oscura pura.
double int_dm (double r, void * params) 
{
  struct param parameters = *(struct param *) params;
  double beta = parameters.beta;  
  double a_dm = parameters.a_dm;
  double int_dm = (1.0 - beta*((R*R)/(r*r)))*((r)/(sqrt(r*r-R*R)))*(((12.0*pow((a_dm+r),4.0)*log((a_dm+r)/(r))) - a_dm*(25.0*a_dm*a_dm*a_dm+52.0*a_dm*a_dm*r+42.0*a_dm*r*r+12.0*r*r*r))/(12.0*pow(a_dm,5)*pow((a_dm+r),4)));
  return int_dm;
}

//Defino la funcion a integrar de la contribucion estelar pura.
double int_ste (double r, void * params) 
{
  struct param parameters = *(struct param *) params;
  double beta = parameters.beta;  
  double a_ste = parameters.a_ste;
  double int_ste = (1.0 - beta*((R*R)/(r*r)))*((r)/(sqrt(r*r-R*R)))*(((12.0*pow((a_ste+r),4.0)*log((a_ste+r)/(r))) - a_ste*(25.0*a_ste*a_ste*a_ste+52.0*a_ste*a_ste*r+42.0*a_ste*r*r+12.0*r*r*r))/(12.0*pow(a_ste,5)*pow((a_ste+r),4)));
  return int_ste;
}

int main(){
  
  FILE *mod,*script;
  int warn;
  double X_s,I_R,s,R,sig_p,rho,re1,re2,result1,error1,result2,error2;
  int Nint = 2000; // Numero de intervalos
  struct param params; // structura de parametros
  gsl_integration_workspace *z = gsl_integration_workspace_alloc(Nint);
  gsl_function F2;
  // parametros del problema (o algo asi)
  //double a = 2.0;
  double a_dm = 2.0;
  double M_dm = 2.0;
  double a_ste = 2.0;
  double M_ste = 2.0;
  double beta = 2.0;

  mod=fopen("model.dat","w");
  //params.a = a;
  params.a_dm = a_dm;
  params.M_dm = M_dm;
  params.a_ste = a_ste;
  params.M_ste = M_ste;
  params.beta = beta;
      
  for(R = 0.01; R < N; R=R+0.01){

    F2.function = &int_dm;
    F2.params = &params;

    gsl_integration_qags(&F2, R, 10000, 0, 1e-4, Nint, z, &result2, &error2); 
    //printf ("result  = % .18f\n", result2);

    if (R < 1.0)
      {      
	
	re2 = result2;
	s = R/a_dm;
	X_s = (1.0 / sqrt(1.0 - s*s)) * log((1.0+(sqrt(1.0-s*s)))/(s));
	I_R = (M_dm/(2.0*M_PI*a_dm*a_dm*GAMMA*(1.0-s*s)*(1.0-s*s)))*((2.0+s*s)*X_s-3.0);
	sig_p = ((2.0*G*M_dm*M_dm*a_dm)/(GAMMA*I_R*2.0*M_PI))*re2;
	//printf("sig_p^2(R) %6lf\n", sig_p);
	fprintf(mod,"%6lf\t %6lf\n", s, sig_p);
      }
    else if(R > 1.01)
      {	       	  
	re2 = result2;
	s = R/a_dm;
	X_s = (1.0 / sqrt(s*s - 1.0)) * acos(1.0/s);
	I_R = (M_dm/(2.0*M_PI*a_dm*a_dm*GAMMA*(1.0-s*s)*(1.0-s*s)))*((2.0+s*s)*X_s-3.0);
	sig_p = ((2.0*G*M_dm*M_dm*a_dm)/(GAMMA*(I_R)*2.0*M_PI))*re2;
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