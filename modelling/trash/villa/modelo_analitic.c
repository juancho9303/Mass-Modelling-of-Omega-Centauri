#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#define N 10
#define M 10e-3
#define GAMMA 1.2
#define G 43007.1

struct param{
  double a;
  double beta;
  double a_ste;
  double a_dm;
  double M_dm;
  double M_ste;
};

double R;

double tot (double r, void * params) 
{
  struct param parameters = *(struct param *) params;
  double a  = parameters.a;
  double a_ste = parameters.a_ste;
  double a_dm = parameters.a_dm;
  double M_ste = parameters.M_ste;
  double M_dm = parameters.M_dm;
  double beta = parameters.beta;
  
  double alfa = ((1.0 - beta*((R*R)/(r*r)))*((r)/(sqrt(r*r-R*R))))/2.0*M_PI;
  
  double M_steM_ste = M_ste*M_ste*a_ste*(((12.0*pow((a_ste+r),4.0)*log((a_ste+r)/(r))) - a_ste*(25.0*a_ste*a_ste*a_ste+52.0*a_ste*a_ste*r+42.0*a_ste*r*r+12.0*r*r*r))/(12.0*pow(a_ste,5)*pow((a_ste+r),4)));
  
/*  double M_steM_dm = M_ste*M_dm*(-(a_ste*(a_ste-a_dm)*a_dm*(a_ste*a_dm*(5.0*a_ste*a_ste-2.0*a_ste*a_dm+5.0*a_dm*a_dm)+4.0*(a_ste+a_dm)*(a_ste*a_ste+a_ste*a_dm+a_dm*a_dm)*r+8.0*((a_ste*a_ste)+(a_ste*a_dm)+(a_dm*a_dm))*r*r+4.0*(a_ste+a_dm)*r*r*r+4.0*(a_ste+r)*(a_ste+r)*(a_dm+r)*(a_dm+r)*(pow((a_ste-a_dm),3.0)*log(r)+a_dm*a_dm*(-3.0*a_ste+a_dm)*log(a_ste+r)-a_ste*a_ste*(a_ste-3.0*a_dm)*log(a_dm+r))))/(2.0*a_ste*a_ste*pow((a_ste-a_dm),3.0)*a_dm*a_dm*(a_ste+r)*(a_ste+r)*(a_dm+r)*(a_dm+r)));
*/

  double M_steM_dm = -M_dm*M_ste*(1.0/(2.0*pow(a_ste,3.0)*pow((a_ste-a_dm),4.0)*a_dm*a_dm*(a_ste+r)*(a_ste+r)*(a_dm+r)))*(2.0*pow((a_ste-a_dm),4)*(a_ste+r)*(a_ste+r)*(a_dm+r)*log(r)-2.0*a_dm*a_dm*(6.0*a_ste*a_ste-4.0*a_ste*a_dm+a_dm*a_dm)*(a_ste+r)*(a_ste+r)*(a_dm+r)*log(a_ste+r)+a_ste*((a_ste-a_dm)*a_dm*(2.0*pow(a_ste,4.0)+4.0*pow(a_ste,3.0)*r-2.0*a_dm*a_dm*r*(a_dm+r)+3.0*a_ste*a_dm*(-a_dm*a_dm+a_dm*r+2.0*r*r)+a_ste*a_ste*(7.0*a_dm*a_dm+7.0*a_dm*r+2.0*r*r))-2.0*a_ste*a_ste*(a_ste-4.0*a_dm)*(a_ste+r)*(a_ste+r)*(a_dm+r)*log(a_dm+r)));
  
  //double M_dmM_ste = ;

  double M_dmM_dm = M_dm*M_dm*a_dm*(((12.0*pow((a_dm+r),4.0)*log((a_dm+r)/(r))) - a_dm*(25.0*a_dm*a_dm*a_dm+52.0*a_dm*a_dm*r+42.0*a_dm*r*r+12.0*r*r*r))/(12.0*pow(a_dm,5)*pow((a_dm+r),4)));
  
  double tot = alfa*(M_steM_ste + M_dmM_dm + M_steM_dm);
  
  return tot;
}

int main(){
  
  FILE *mod,*script;
  int warn;
  double X_s,I_R,s,R,sig_p,rho,re1,re2,result1,error1,result2,error2;
  int Nint = 1000; // Numero de intervalos
  struct param params; // structura de parametros
  gsl_integration_workspace *z = gsl_integration_workspace_alloc(Nint);
  gsl_function F1;
  // parametros del problema (o algo asi)
  double a = 1.0;
  double a_dm = 1.0;
  double a_ste = 1.0;
  double M_ste = 2.0;
  double M_dm = 2.0;
  double beta = 2.0;

  mod=fopen("model.dat","w");
  params.a = a;
  params.beta = beta;
  params.a_ste = a_ste;
  params.a_dm = a_dm;
  params.M_ste = M_ste;
  params.M_dm = M_dm;
      
  for(R = 0.01; R < N; R=R+0.01){

    F1.function = &tot;
    F1.params = &params;
    
    // Para todos es la misma integral
    gsl_integration_qags(&F1, R, 1000, 0, 1e-4, Nint, z, &result1, &error1); 
    printf ("result  = % .18f\n", result1);

    if (R < 1.0)
      {      
	
	re1 = result1;
	s = R/a;
	X_s = (1.0 / sqrt(1.0 - s*s)) * log((1.0+(sqrt(1.0-s*s)))/(s));
	I_R = (M/(2.0*M_PI*a*a*GAMMA*(1.0-s*s)*(1.0-s*s)))*((2.0+s*s)*X_s-3.0);
	sig_p = ((2.0*G*M*M*a)/(GAMMA*I_R*2.0*M_PI))*re1;
	//printf("sig_p^2(R) %6lf\n", sig_p);
	fprintf(mod,"%6lf\t %6lf\n", s, sig_p);
      }
    else if(R > 1.01)
      {	       	  
	re1 = result1;
	s = R/a;
	X_s = (1.0 / sqrt(s*s - 1.0)) * acos(1.0/s);
	I_R = (M/(2.0*M_PI*a*a*GAMMA*(1.0-s*s)*(1.0-s*s)))*((2.0+s*s)*X_s-3.0);
	sig_p = ((2.0*G*M*M*a)/(GAMMA*(I_R)*2.0*M_PI))*re1;
	//printf("sig_p^2(R) %6lf\n", sig_p);
	fprintf(mod,"%6lf\t %6lf\n", s, sig_p);
      }      
  }

  gsl_integration_workspace_free(z);
  
  fclose(mod);
  
  script = fopen( "script.gpl", "w" );
  fprintf(script, "set terminal png\nset output 'sigma.png'\nset nokey\n");
  fprintf( script, "set title 'Sigma Proyectada vs R'\n" );
  fprintf( script, "plot 'model.dat' u 1:2 w d\n");
  fclose(script);
  warn = system("gnuplot script.gpl");
  
  return(warn);
}
