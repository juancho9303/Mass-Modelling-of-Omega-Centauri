#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#define N 10
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

double R,result1,result2,result3,result4,error1,error2,error3,error4,resultt,errort;
int Nint = 2000; // Numero de intervalos

//M*M*
double M_steM_ste (double r, void * params) 
{
  struct param parameters = *(struct param *) params;
  double M_ste = parameters.M_ste;
  double a_ste = parameters.a_ste;
  double M_steM_ste = (M_ste*M_ste*a_ste)/(r*pow((r+a_ste),5.0));
  return M_steM_ste;
}

//M*M_dm
double M_steM_dm (double r, void * params) 
{
  struct param parameters = *(struct param *) params;
  double M_ste = parameters.M_ste;
  double a_ste = parameters.a_ste;
  double M_dm = parameters.M_dm;
  double a_dm = parameters.a_dm;
  double M_steM_dm = (M_ste*M_dm*a_ste)/(r*pow((r+a_ste),3.0)*pow((r+a_dm),2.0));
  return M_steM_dm;
}

//M_dmM*
double M_dmM_ste (double r, void * params) 
{
  struct param parameters = *(struct param *) params;
  double M_ste = parameters.M_ste;
  double a_ste = parameters.a_ste;
  double M_dm = parameters.M_dm;
  double a_dm = parameters.a_dm;
  double M_dmM_ste = (M_dm*M_ste*a_dm)/(r*pow((r+a_dm),3.0)*pow((r+a_ste),2.0));
  return M_dmM_ste;
}

//M_dmM_dm
double M_dmM_dm (double r, void * params) 
{
  struct param parameters = *(struct param *) params;
  double M_dm = parameters.M_dm;
  double a_dm = parameters.a_dm;
  double M_dmM_dm = (M_dm*M_dm*a_dm)/(r*pow((r+a_dm),5.0));
  return M_dmM_dm;
}

//Integral total
double int_tot (double r, void * params) 
{
  struct param parameters = *(struct param *) params;
  double beta = parameters.beta;  
  double a_dm = parameters.a_dm;
  double M_dm = parameters.a_dm;
  double a_ste = parameters.a_dm;
  double M_ste = parameters.a_dm;
  //Cuando hay dark matter y stellar
  
  // esto deberia venir aca
  
  gsl_integration_workspace *j = gsl_integration_workspace_alloc(Nint);
  gsl_function F1;
  gsl_function F2;
  gsl_function F3;
  gsl_function F4;  
  
  F1.function = &M_steM_ste;
  F1.params = &params;
  gsl_integration_qags(&F1, r, 600, 0, 1e-4, Nint, j, &result1, &error1); 
  printf ("result  = % .18f\n", result1);
  
  F2.function = &M_steM_dm;
  F2.params = &params;
  gsl_integration_qags(&F2, r, 600, 0, 1e-4, Nint, j, &result2, &error2); 
  //printf ("result  = % .18f\n", result2);
  
  F3.function = &M_dmM_ste;
  F3.params = &params;
  gsl_integration_qags(&F3, r, 600, 0, 1e-4, Nint, j, &result3, &error3); 
  //printf ("result  = % .18f\n", result3);
  
  F4.function = &M_dmM_dm;
  F4.params = &params;
  gsl_integration_qags(&F4, r, 600, 0, 1e-4, Nint, j, &result4, &error4); 
  //printf ("result  = % .18f\n", result4);
    
  double int_tot = (1.0 - beta*((R*R)/(r*r)))*((r)/(sqrt(r*r-R*R)))*(result1+result2+result3+result4);
  return int_tot;
} 

int main(){
   
  FILE *mod,*script;
  int warn;
  double X_s,I_R,s,R,sig_p,rho,re1,re_tot;
  //int Nint = 2000; // Numero de intervalos
  struct param params; // structura de parametros
  
  //printf ("int_tot = % .18f\n", int_tot);
  
  /* gsl_integration_workspace *j = gsl_integration_workspace_alloc(Nint);
     gsl_function F1;
     gsl_function F2;
     gsl_function F3;
     gsl_function F4; */
  
  gsl_integration_workspace *z = gsl_integration_workspace_alloc(Nint);
  gsl_function Ft;
  
  double a_dm = 1.0;
  double M_dm = 2.0;
  double a_ste = 1.0;
  double M_ste = 2.0;
  double beta = 2.0;
  
  mod=fopen("model.dat","w");
  

    params.a_dm = a_dm;
    params.M_dm = M_dm;
    params.a_ste = a_ste;
    params.M_ste = M_ste;
    params.beta = beta;
 /*   
    F1.function = &M_steM_ste;
    F1.params = &params;
    gsl_integration_qags(&F1, 0.001, 600, 0, 1e-4, Nint, j, &result1, &error1); 
    printf ("result  = % .18f\n", result1);
    
    F2.function = &M_steM_dm;
    F2.params = &params;
    gsl_integration_qags(&F2, 0.001, 600, 0, 1e-4, Nint, j, &result2, &error2); 
    printf ("result  = % .18f\n", result2);
    
    F3.function = &M_dmM_ste;
    F3.params = &params;
    gsl_integration_qags(&F3, 0.001, 600, 0, 1e-4, Nint, j, &result3, &error3); 
    printf ("result  = % .18f\n", result3);
    
    F4.function = &M_dmM_dm;
    F4.params = &params;
    gsl_integration_qags(&F4, 0.001, 600, 0, 1e-4, Nint, j, &result4, &error4); 
    printf ("result  = % .18f\n", result4);
  */
 
  /*   //Integral total
       double int_tot (double r, void * params) 
       {
       struct param parameters = *(struct param *) params;
       double beta = parameters.beta;  
       double a_dm = parameters.a_dm;
       //double int_tot = (1.0 - beta*((R*R)/(r*r)))*((r)/(sqrt(r*r-R*R)))*(result1+result2+result3+result4);
       double int_tot = (1.0 - beta*((R*R)/(r*r)))*((r)/(sqrt(r*r-R*R)))*(((12.0*pow((a_dm+r),4.0)*log((a_dm+r)/(r))) - a_dm*(25.0*a_dm*a_dm*a_dm+52.0*a_dm*a_dm*r+42.0*a_dm*r*r+12.0*r*r*r))/(12.0*pow(a_dm,5)*pow((a_dm+r),4)));
       return int_tot;
       } */  
  
  for(R = 0.01; R < N; R=R+0.01){
    
    Ft.function = &int_tot;
    Ft.params = &params;
    gsl_integration_qags(&Ft, R, 1000, 0, 1e-4, Nint, z, &resultt, &errort); 
    //printf ("result  = % .18f\n", resultt);
    
    if (R < 1.0)
      {      	
	re_tot = resultt;
	s = R/a_ste;
	X_s = (1.0 / sqrt(1.0 - s*s)) * log((1.0+(sqrt(1.0-s*s)))/(s));
	I_R = (M_ste/(2.0*M_PI*a_ste*a_ste*GAMMA*(1.0-s*s)*(1.0-s*s)))*((2.0+s*s)*X_s-3.0);
	sig_p = ((2.0*G)/(GAMMA*I_R*2.0*M_PI))*re_tot;
	//printf("sig_p^2(R) %6lf\n", sig_p);
	fprintf(mod,"%6lf\t %6lf\n", R, sig_p);
      }
    else if(R > 1.01)
      {	       	  
	re_tot = resultt;
	s = R/a_dm;
	X_s = (1.0 / sqrt(s*s - 1.0)) * acos(1.0/s);
	I_R = (M_dm/(2.0*M_PI*a_dm*a_dm*GAMMA*(1.0-s*s)*(1.0-s*s)))*((2.0+s*s)*X_s-3.0);
	sig_p = ((2.0*G)/(GAMMA*(I_R)*2.0*M_PI))*re_tot;
	//printf("sig_p^2(R) %6lf\n", sig_p);
	fprintf(mod,"%6lf\t %6lf\n", R, sig_p);
      }      
  }
  
  // gsl_integration_workspace_free(j);
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
