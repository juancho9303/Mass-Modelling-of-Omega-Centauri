#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

/*
  SE UTILIZA UN MODELO DE HERNQUIST MODIFICADO PARA ENCONTRAR LA
  FORMA DE LA DISPERSION DE VELOCIDADES PROYECTADA Y PODER AJUSTAR
  LOS PARAMETROS QUE MEJOR SE AJUSTEN A LAS OBSERVACIONES   
*/

#define K 90           //90 parsecs que esta por encima de los ~60pc de omega centauri
//#define GAMMA 1.0        //razon masa luminosidad
#define G 43007.1        //unidades de Gadget
#define EPSILON 1     //para suavizar denominadores para evitar indeterminaciones

/* UNIDADES:
   
   UnitMass_in_g = 1.989e+43 
   UnitTime_in_s = 3.08568e+16 
   UnitVelocity_in_cm_per_s = 100000 
   Unidad de Masa = 10^10 Msun/h
   Unidad de longitud = 1kpc/h
*/

// GENERO LA ESTRUCTURA DE PARAMETROS QUE VOY A AJUSTAR
struct param
{
  double beta, a, M;  
};

double R;

// ACA SE DEFINE LA INTEGRAL QUE CONTIENE TODOS LOS TERMINOS DE LAS CONTRIBUCIONES DE MASA ESTELAR Y MATERIA OSCURA Y QUE DEPENDE DE r 

double integrando(double r, void * params) 
{
  struct param parameters = *(struct param *) params;
  double a  = parameters.a;
  double beta = parameters.beta;  
  double integrando = (1.0 - beta*((R*R)/(r*r)))*((r)/(sqrt(r*r-R*R)))*(((12.0*pow((a+r),4.0)*log((a+r)/(r))) - a*(25.0*a*a*a+52.0*a*a*r+42.0*a*r*r+12.0*r*r*r)))/(12.0*pow(a,5)*pow((a+r),4));
  return integrando;
}

double integrando2(double r, void * params) 
{
  struct param parameters = *(struct param *) params;
  double a_min  = parameters.a;
  double beta = parameters.beta;  
  double integrando2 = (1.0 - beta*((R*R)/(r*r)))*((r)/(sqrt(r*r-R*R)))*(((12.0*pow((a_min+r),4.0)*log((a_min+r)/(r))) - a_min*(25.0*a_min*a_min*a_min+52.0*a_min*a_min*r+42.0*a_min*r*r+12.0*r*r*r)))/(12.0*pow(a_min,5)*pow((a_min+r),4));
  return integrando2;
}

// AHORA INICIO EL PROGRAMA PRINCIPAL QUE VA A HACER TODOS LOS CALCULOS PARA ENCONTRAR LA DSPERSION DE VELOCIDADES PROYECTADA

int main()
{
  
  // DEFINO LOS PUNTEROS DE LOS ARCHIVOS, LAS VARIABLES DE TIPO ENTERO Y DOUBLE Y LA ALOCACION DE MEMORIA DE LA RUTINA DE GSL QUE VOY A USAR PARA LA INTEGRAL
  
  FILE *mod,*script,*sigma,*chi_cua;
  int warn, t, i, j, k;
  int Nint = 1000; // Numero de intervalos
  double X_s, I_R, s, R[90], R_arcmin[90], sig_p[90], rho, re1, re2, result1, error1, a, M, beta, R_o[12], sig_p_o[12], chi2, chi_t, a_min, M_min, R_arcmin_i, sig_p_i;
  double chi_min = 1000000.0;
  double gamma, gamma_min;
  struct param params;
  gsl_integration_workspace *z = gsl_integration_workspace_alloc(Nint);
  gsl_function F1;
  gsl_function F2;
  gsl_interp_accel *acc  =   gsl_interp_accel_alloc ();
  gsl_spline *spline     =   gsl_spline_alloc (gsl_interp_cspline, 90);
  
  sigma=fopen("sigma.dat","r"); 
  
  // ESTE CICLO ES PARA LEER EL ARCHIVO DE DATOS CON LOS DATOS OBSERVACIONALES DE SIGMA
  
  for(j=0.0; j<12; j++)
    {        
      fscanf(sigma, "%lf %lf", &R_o[j], &sig_p_o[j]);
    }      
  
  fclose(sigma);
  
  // HAGO LAS VARIACIONES DE LOS PARAMETROS
  
  beta  =  0.0001;
  //  a = 0.055;
  //  M = 0.0005;
  //  gamma = 1.5;
  
  //mod=fopen("model.dat","w"); 
  chi_cua=fopen("chi_cuadrado.dat","w");  
  
  for(gamma = 1.5; gamma <= 3.0;  gamma = gamma + 0.5)
    { 
      for(a = 0.005; a <= 0.06;  a = a + 0.01)
	{
	  for(M = 0.00001; M <= 0.0009;  M = M + 0.0001)
	    { 
	      params.beta  =  beta;
	      params.a   =  a;
	      params.M   =  M;
	      
	      // ESTE CICLO ES PARA R (RADIO PROYECTO) YA QUE SIGMA DEPENDE DE R, ADEMAS VA A SER EL LIMITE INFERIOR DE LA INTEGRAL
	      
	      for(i = 0; i < K; i ++)          
		{
		  R[i] = (i+1.0)/(1000.0);
		  F1.function = &integrando;
		  F1.params = &params;
		  //printf("%lf\n",integrando);
		  gsl_integration_qags(&F1, R[i], 0.1, 0, 1e-4, Nint, z, &result1, &error1); 
		  
		  re1 = result1;
		  printf("%d %lf\n", i, re1);
		  s = R[i]/a;
		  
		  // ES NECESARIO PONER ESTOS CONDICIONALES YA QUE LA FUNCION X(S) ESTA DEFINIDA A TRAMOS
		  
		  if (s < 1.0 - 1.0e-9)
		    {      	
		      X_s = (1.0 / sqrt(1.0 - s*s)) * log((1.0+(sqrt(1.0-s*s)))/(s));
		      I_R = (M/(2.0*M_PI*a*a*gamma*(1.0-s*s)*(1.0-s*s)))*((2.0+s*s)*X_s-3.0);
		    }
		  
		  if(s >= 1.0-1.0e-10 && s <= 1.0+1.0e-10)
		    if(s = 1.0)
		      {
			X_s = 1.0;
			I_R = (2.0*M)/(15.0*M_PI*a*a*gamma);
		      }
		  
		  if (s > 1.0 + 1.0e-9)
		    {	       	  
		      X_s = (1.0 / sqrt(s*s - 1.0)) * acos(1.0/s);
		      I_R = (M/(2.0*M_PI*a*a*gamma*(1.0-s*s)*(1.0-s*s)))*((2.0+s*s)*X_s-3.0);
		    }      
		  
		  // SE HACEN LOS CALCULOS UNA VEZ SE TIENE EL VALOR DE LA INTEGRAL
		  
		  sig_p[i] = sqrt(((M*M*G*a)/(gamma*(I_R)*M_PI))*re1);
		  //CONVIERTO DE KILOPARSECS A MINUTOS DE ARCO
		  R_arcmin[i] = R[i]/(4.80839)*3437.75;
		  //fprintf (chi_cua, "%g %g\n", R_arcmin[i], sig_p[i]);
		}
	      fprintf(chi_cua,"\n");
	      
	      gsl_spline_init (spline, R_arcmin, sig_p, 90);
	      
	      for(R_arcmin_i = R_arcmin[0]; R_arcmin_i < R_arcmin[89]; R_arcmin_i += 0.001)
		{
		  sig_p_i = gsl_spline_eval (spline, R_arcmin_i, acc);
		  //printf("%lf %lf\n",R_arcmin_i, sig_p_i);
		  //printf("%lf\n",sig_p_i);
		  
		  for(j=0.0; j<12.0; j=j+1.0)
		    {
		      //printf ("%g %g\n", R_o[j], sig_p_o[j]);
		      if( R_arcmin_i - R_o[j] < 0.001 && R_arcmin_i - R_o[j] > -0.001 )
			//if( R_arcmin_i = R_o[j] )
			{
			  if(sig_p_i > 0.0)
			    {
			      //printf("%lf %lf\n",sig_p_i,sig_p_o[j]);
			      chi2 = (sig_p_i-sig_p_o[j])*(sig_p_i-sig_p_o[j]);
			      //printf("%lf\n",chi2);
			      chi_t = chi_t + chi2;
			      //printf("%lf\n",chi_t);
			    }       
			}
		    }
		}
	      
	      if(chi_t < chi_min)
		{
		  chi_min   =  chi_t;
		  a_min     =  a;
		  M_min     =  M;
		  gamma_min =  gamma;
		  //printf("%lf %lf\n", a_min, M_min);
		  
		  params.beta  =  beta;
		  params.a     =  a_min;
		  params.M     =  M_min;  
		  
		}
	      //printf("%lf\n", chi_t);
	      chi_t = 0.0;      
	    }
	}
    }
  
  params.a   =  a_min;
  params.M   =  M_min;
  
  printf("%g %g\n", a_min, M_min);
  
  for(i = 0; i < K; i ++)          
    { 
      R[i] = (i+1.0)/(1000.0);
      F2.function = &integrando2;
      F2.params = &params;
      gsl_integration_qags(&F2, R[i], 0.1, 0, 1e-4, Nint, z, &result1, &error1); 
      re1 = result1;
      s = R[i]/a_min;
      
      if (s < 1.0 - 1.0e-9)
	{      	
	  X_s = (1.0 / sqrt(1.0 - s*s)) * log((1.0+(sqrt(1.0-s*s)))/(s));
	  I_R = (M_min/(2.0*M_PI*a_min*a_min*gamma_min*(1.0-s*s)*(1.0-s*s)))*((2.0+s*s)*X_s-3.0);
	}
      
      if(s >= 1.0-1.0e-10 && s <= 1.0+1.0e-10)
	if(s = 1.0)
	  {
	    X_s = 1.0;
	    I_R = (2.0*M_min)/(15.0*M_PI*a_min*a_min*gamma_min);
	  }
      
      if (s > 1.0 + 1.0e-9)
	{	       	  
	  X_s = (1.0 / sqrt(s*s - 1.0)) * acos(1.0/s);
	  I_R = (M_min/(2.0*M_PI*a_min*a_min*gamma_min*(1.0-s*s)*(1.0-s*s)))*((2.0+s*s)*X_s-3.0);
	}      
      
      sig_p_i = sqrt(((G*M_min*M_min*a_min)/(gamma_min*(I_R)*M_PI))*re1);
      R_arcmin_i = R[i]/(4.80839)*3437.75;
      fprintf(chi_cua,"%16.8e\t %16.8e\n", R_arcmin_i, sig_p_i);
    }
  
  // MUESTRO CUALES SON LOS PARAMETROS QUE MEJOR SE AJUSTAN A LOS DATOS OBSERVACIONALES
  
  printf("xi=%lf\t a=%lf\t M=%lf gamma=%lf\n",chi_min,a_min,M_min,gamma_min);
  //fclose(mod);

fclose(chi_cua);
gsl_integration_workspace_free(z);
gsl_spline_free (spline);
gsl_interp_accel_free (acc);

script = fopen( "script.gpl", "w" );
fprintf( script, "plot 'sigma.dat'pt 7 ps 1.2\n" );
fprintf(script, "set grid\nset terminal png\nset output 'sigma.png'\nset nokey\n");
fprintf( script, "set title 'Sigma Proyectada vs R'\n" );
//fprintf( script, "set xrange [0:40]\n" );
//fprintf( script, "set yrange [4:14]\n" );
fprintf( script, "set xlabel 'Radius in Arcmin'\n" );
fprintf( script, "set ylabel 'Projected velocity dispersion in km/s'\n" );
fprintf( script, "replot 'chi_cuadrado.dat' u 1:2 w l\n");
fclose(script);

warn = system("gnuplot script.gpl");

return(warn);
}
