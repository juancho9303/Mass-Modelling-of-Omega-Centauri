#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

/*
  SE UTILIZA UN MODELO DE HERNQUIST MODIFICADO PARA ENCONTRAR LA FORMA DE LA DISPERSION DE VELOCIDADES PROYECTADA Y PODER AJUSTAR LOS PARAMETROS QUE MEJOR SE AJUSTEN A LAS OBSERVACIONES   
*/

#define K 90                 //90 parsecs que esta por encima de los ~60pc de omega centauri
#define G 43007.1            //unidades de Gadget
#define EPSILON 1            //para suavizar denominadores para evitar indeterminaciones
#define DISTANCE 4.80839
#define LIMIT_RADIUS 0.5
#define ZERO 1.9e-10

// GENERO LA ESTRUCTURA DE PARAMETROS QUE VOY A AJUSTAR
struct param
{
  double beta, a_s, a_dm, M_s, M_dm, R;  
};

#include "routines.c"

int evaluate_integral(double R, double beta, double gamma, double a_dm, double a_s, double M_dm, double M_s, double *sig_p, double *R_arcmin)
{
  
  int Nint = 1000; // Numero de intervalos
  double result, error, s, aux, X_s, I_R;
  
  gsl_integration_workspace *z = gsl_integration_workspace_alloc(Nint);
  gsl_function F1;
  
  struct param params;
  
  params.beta  =  beta;
  params.a_s   =  a_s;
  params.a_dm  =  a_dm;
  params.M_s   =  M_s;
  params.M_dm  =  M_dm;       
  params.R     =  R;
  
  F1.function = &integrando;
  F1.params = &params;
  
  gsl_integration_qags(&F1, R, LIMIT_RADIUS, 0, 1e-6, Nint, z, &result, &error); 
  
  s = R/a_s;
  
  if (s < (1.0-ZERO))
    { 
      aux = 1.0 - s*s;
      X_s = (1.0 / sqrt(aux)) * log((1.0 + sqrt(aux))/s);
      I_R = (M_s/(2.0*M_PI*a_s*a_s*gamma*aux*aux))*( (2.0+s*s)*X_s - 3.0 );
    }
  
  if(s >= (1.0-ZERO) && s <= (1.0+ZERO) )
    {
      X_s = 1.0;
      I_R = (2.0*M_s)/(15.0*M_PI*a_s*a_s*gamma);
    }
  
  if (s > (1.0 + ZERO))
    {
      X_s = (1.0 / sqrt(s*s - 1)) * acos(1.0/s);
      I_R = (M_s/(2.0*M_PI*a_s*a_s*gamma*(1.0-s*s)*(1.0-s*s)))*((2.0+s*s)*X_s - 3.0);
    }      
  
  // SE HACEN LOS CALCULOS UNA VEZ SE TIENE EL VALOR DE LA INTEGRAL
  
  *sig_p = sqrt( (G*result) / (gamma*I_R*M_PI) ) ;
  
  //CONVIERTO DE KILOPARSECS A MINUTOS DE ARCO
  *R_arcmin = R*3437.75/DISTANCE;
  
  gsl_integration_workspace_free(z);
  
  return 0;  
}

// AHORA INICIO EL PROGRAMA PRINCIPAL QUE VA A HACER TODOS LOS CALCULOS PARA ENCONTRAR LA DSPERSION DE VELOCIDADES PROYECTADA

int main(int argc, char *argv[])
{
  
  int warn, i, j, nread, NUMBER_OF_SIGMA_BINS;
  
  double X_s, I_R, s, R[90], R_arcmin[90], sig_p[90], a_dm, a_s, M_s, M_dm, beta;
  double R_o[12], sig_p_o[12], numb_error[12], chi2, a_dm_min, a_s_min, M_dm_min, M_s_min, beta_min, sig_p_i;
  double chi_min = 1000000.0;
  double gamma, gamma_min, aux;
  char *infile;
  FILE *script,*sigma,*chi_cua;
  
  double gamma_ini, gamma_end, Delta_gamma;
  double adm_ini, adm_end, Delta_adm;
  double as_ini, as_end, Delta_as;
  double Mdm_ini, Mdm_end, Delta_Mdm;
  double Ms_ini, Ms_end, Delta_Ms;
  double beta_ini, beta_end, Delta_beta;
  
  
  infile = argv[1];                          // input file with velocity dispersion bins!
  NUMBER_OF_SIGMA_BINS = atoi(argv[2]);      // eventually 12!
  
  if(argc != 3)
    {
      printf(" Execute as ./exec <infile> <number of lines>\n");
      exit(0);
    }
  
  sigma=fopen(infile,"r"); 
  for(j=0.0; j<NUMBER_OF_SIGMA_BINS; j++)
    nread = fscanf(sigma,"%lf %lf %lf", &R_o[j], &sig_p_o[j], &numb_error[j]);
  fclose(sigma);
  
  a_s =  0.00223;  //Stellar scalength of 2.234 parsecs from the fitting of Noyola's data
  
  gamma_ini  = 0.1;        gamma_end  = 3.0;      Delta_gamma  = 0.5;
  adm_ini    = 0.001;      adm_end    = 0.06;     Delta_adm    = 0.01;
  Mdm_ini    = 0.00001;    Mdm_end    = 0.0008;   Delta_Mdm    = 0.0001;
  Ms_ini     = 0.00001;    Ms_end     = 0.0008;   Delta_Ms     = 0.0001;
  beta_ini   = 0.00001;    beta_end   = 1.0;      Delta_beta   = 0.03;
  
  // HAGO LAS VARIACIONES DE LOS PARAMETROS
  
  int Ntotal_iterations, counter, cuenta;
  Ntotal_iterations = (int) ( ((beta_end-beta_ini)/Delta_beta)*((gamma_end-gamma_ini)/Delta_gamma)*((adm_end-adm_ini)/Delta_adm)*((Mdm_end-Mdm_ini)/Delta_Mdm)*((Ms_end-Ms_ini)/Delta_Ms) );
  
  printf("  Running %d iterations in parameter space\n", Ntotal_iterations);
  
  counter=0;
  
  /*DEFINICION DEL VECTOR R[i]*/
  double ii = 0.0;
  for(i=0;i<K;i++)
    {
      R[i] = (ii+1.0)/1000.0;
      if(i==0)
	R[i] = R[i]*0.5;
      ii = ii +1.0;
    }
  
  for(beta = beta_ini; beta <= beta_end; beta = beta + Delta_beta)
    { 
      for(gamma = gamma_ini; gamma <= gamma_end;  gamma = gamma + Delta_gamma)
	{     
	  for(a_dm = adm_ini; a_dm <= adm_end;  a_dm = a_dm + Delta_adm)
	    {
	        // PARA EVITAR INDETERMINACIONES EN ALGUNOS DENOMINADORES
		  if (a_s != a_dm)
		    {
		      for(M_dm = Mdm_ini; M_dm < Mdm_end;  M_dm = M_dm + Delta_Mdm)
			{
			  for(M_s = Ms_ini; M_s < Ms_end;  M_s = M_s + Delta_Ms)
			    {
			      // ESTE CICLO ES PARA R (RADIO PROYECTADO) YA QUE SIGMA DEPENDE DE R, ADEMAS VA A SER EL LIMITE INFERIOR DE LA INTEGRAL
			      
			      //printf("%16.8lf\n",a_dm);
			      
			      for(i = 0; i < K; i ++)          
				{
				  evaluate_integral(R[i], beta, gamma, a_dm, a_s, M_dm, M_s, &sig_p[i], &R_arcmin[i]);
				  
				}

			      gsl_interp_accel *acc  =   gsl_interp_accel_alloc ();
			      gsl_spline *spline     =   gsl_spline_alloc (gsl_interp_cspline, 90);

			      gsl_spline_init (spline, R_arcmin, sig_p, 90);
			      
			      chi2 = 0.0;
			      for(j=0; j<NUMBER_OF_SIGMA_BINS; j++)
				{
				  sig_p_i = gsl_spline_eval (spline, R_o[j], acc);
				  chi2 += ( (sig_p_i-sig_p_o[j])*(sig_p_i-sig_p_o[j]) ) / numb_error[j];
				}

			      gsl_spline_free (spline);
			      gsl_interp_accel_free (acc);
			      
			      if(chi2 < chi_min)
				{
				  chi_min   = chi2;
				  a_dm_min  = a_dm;
				  a_s_min   = a_s;
				  M_dm_min  = M_dm;
				  M_s_min   = M_s;
				  gamma_min = gamma;
				  beta_min = beta;
				}
			      
			      if((counter%1000) == 0)
				printf("Ready %d iterations of %d\n",counter, Ntotal_iterations);
			      counter++;
			      
			    }// for Ms
			}// for Mdm	      
		    }// if for the singularity
	    }// for adm
	}//for gamma
    }//for beta
  
  //printf ("\n");
  
  chi_cua=fopen("chi_cuadrado.dat","w"); 
  
  for(i = 0; i < K; i ++)          
    {
      
      evaluate_integral(R[i], beta_min, gamma_min, a_dm_min, a_s_min, M_dm_min, M_s_min, &sig_p[i], &R_arcmin[i]);
      fprintf(chi_cua,"%16.8e %16.8e\n", R_arcmin[i], sig_p[i]);
      
    }
  
  fclose(chi_cua);
  
  // MUESTRO CUALES SON LOS PARAMETROS QUE MEJOR SE AJUSTAN A LOS DATOS OBSERVACIONALES
  
  printf("xi=%16.8lf beta=%16.8lf a_dm/pc=%16.8lf a_s/pc=%16.8lf M_dm/Msun=%16.8elf M_s/Msun=%16.8e gamma=%lf\n",
	 chi_min, beta_min, 1000*a_dm_min, 1000*a_s_min, (1.0e10)*M_dm_min, (1.0e10)*M_s_min, gamma_min);
  
  
  script = fopen( "script.gpl", "w" );
  fprintf( script, "plot 'sigma.dat'pt 7 ps 1.2\n" );
  fprintf( script, "replot 'chi_cuadrado.dat' u 1:2 w l\n");
  fprintf( script, "reset\n");
  fprintf(script, "set grid\nset terminal png\nset output 'sigma.png'\nset nokey\n");
  fprintf( script, "set title 'Sigma Proyectada vs R'\n" );
  fprintf( script, "set xrange [0:40]\n" );
  fprintf( script, "set yrange [4:14]\n" );
  fprintf( script, "set xlabel 'Radius in Arcmin'\n" );
  fprintf( script, "set ylabel 'Projected velocity dispersion in km/s'\n" );
  fprintf( script, "plot 'sigma.dat'pt 7 ps 1.2\n" );
  fprintf( script, "replot 'chi_cuadrado.dat' u 1:2 w l\n");
  fprintf( script, "pause -1\n");
  fclose(script);
  
  warn = system("gnuplot script.gpl");
  
  return(warn);
}
