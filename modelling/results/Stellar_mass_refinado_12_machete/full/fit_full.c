#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

//////////////////////////////////////////////////////////////////////////////
// HERNQUIST MODEL FOR THE PROJECTED VELOCITY DISPERSION IN OMEGA CENTAURI  //
//////////////////////////////////////////////////////////////////////////////

/*  SISTEM OF UNITS 

  1 unit of mass   = 1.0e5Msun
  1 unit of lenght = 1.0 pc
  1 unit velocity  = 1km/s
  G = 430.071
  
  Hubble (internal units) = 0.0001
  UnitMass_in_g = 1.989e+38 
  UnitTime_in_s = 3.08568e+13 
  UnitVelocity_in_cm_per_s = 100000 
  UnitDensity_in_cgs = 6.76991e-18 
  UnitEnergy_in_cgs = 1.989e+48 
  
*/

#define K 90                 //90 PARSECS, OUT OF THE ~60 PC OF OMEGA CENTAURI
#define G 430.071            
#define EPSILON 1            // TO SOFTEN INDETERMINATIONS
#define DISTANCE 4808.39
#define LIMIT_RADIUS 1500
#define ZERO 1.9e-10

struct param
{
  double beta, a_s, a_dm, M_s, M_dm, R;  
};

// THE FILE ROUTINES.C HAS THE LARGE EQUATIONS OF THE INTEGRAND 

#include "routines.c"

//THIS IS THE FUNCTION THAT CONTAINS ALL THE CALCULATIONS, INCLUDING THE INTEGRALS AND INTERPOLATIONS 

int evaluate_integral(double R, double beta, double gamma, double a_dm, double a_s, double M_dm, double M_s, double *sig_p, double *R_arcmin)
{
  
  int Nint = 10000;    // NUMBER OF INTERVALS
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

  //////////////////////////////////////////////////////////////
  //  THIS IS JUST TO PRINT THE INTEGRAND TO SEE ITS BEHAVIOR  //
  ///////////////////////////////////////////////////////////////
  /* 
    {
    FILE *pf=NULL;
    double r;
    pf=fopen("puto_archivo.dat","w");
    for(r=0; r<LIMIT_RADIUS; r=r+1)
    {
    fprintf(pf,"%16.8e %16.8e\n",r, integrando(r, &params));
    }
    fclose(pf);
    exit(0);
    } 
  */
  ///////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  
  F1.function = &integrando;
  F1.params = &params;
  
  //gsl_integration_qags(&F1, R, LIMIT_RADIUS, 1.0e-8, 1e-6, Nint, z, &result, &error); 
  gsl_integration_qag(&F1, R, LIMIT_RADIUS, 1.0e-5, 1e-5, Nint, 1, z, &result, &error); 
  //printf("R = %lf I=%lf\n",R, result);
  
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
  
  *sig_p = sqrt( (G*result) / (gamma*I_R*M_PI) ) ;
  
  //THIS LINE CONVERTS FROM PARSECS TO ARCMIN SINCE THE OBSERVATIONAL DATA IS IN THOSE UNITS
  *R_arcmin = R*3437.75/DISTANCE;
  
  gsl_integration_workspace_free(z);
  
  return 0;  
}

// THE MAIN PROGRAM WILL RUN FOR ALL THE COMBINATIONS OF PARAMETERS AND FIT THE XI SQUARE TO GET THE BEST COMBINATION

int main(int argc, char *argv[])
{
  
  int warn, i, j, nread, NUMBER_OF_SIGMA_BINS; 
  double X_s, I_R, s, R[90], R_arcmin[90], sig_p[90], a_dm, a_s, M_s, M_dm, beta;
  double R_o[12], sig_p_o[12], numb_error[12], chi2, a_dm_min, a_s_min, M_dm_min, M_s_min, beta_min, sig_p_i;
  double chi_min = 100000000000000000.0;
  double gamma, gamma_min, aux;
  char *infile;
  FILE *script,*sigma,*best_model;
  
  double gamma_ini, gamma_end, Delta_gamma;
  double adm_ini, adm_end, Delta_adm;
  double as_ini, as_end, Delta_as;
  double Mdm_ini, Mdm_end, Delta_Mdm;
  double Ms_ini, Ms_end, Delta_Ms;
  double beta_ini, beta_end, Delta_beta;
  
  infile = argv[1];                          // INPUT FILE WITH VELOCITY DISPERSION BINS
  NUMBER_OF_SIGMA_BINS = atoi(argv[2]);      // EVENTUALLY 12!
  
  if(argc != 3)
    {
      printf(" Execute as ./exec <infile> <number of lines>\n");    // THIS WILL TELL HOW TO PROPERLY DO THE EXECUTION
      exit(0);
    }
  
  // READ THE FILE WITH THE OBSERVATIONAL DATA
  sigma=fopen(infile,"r");                      
  for(j=0.0; j<NUMBER_OF_SIGMA_BINS; j++)
    nread = fscanf(sigma,"%lf %lf %lf", &R_o[j], &sig_p_o[j], &numb_error[j]);
  fclose(sigma);
  
  M_s = 49.6;

  // THE VARIATIONS OF THE PARAMETERS

  gamma_ini  = 0.6;        gamma_end  = 0.8;      Delta_gamma  = 0.06;
  adm_ini    = 12.0;        adm_end    = 14.0;     Delta_adm    = 0.06;
  as_ini     = 30.0;        as_end     = 32.0;     Delta_as     = 0.06;
  Mdm_ini    = 18.0;        Mdm_end    = 20.0;     Delta_Mdm    = 0.06;
  beta_ini   = 0.8;     beta_end   = 1.0;      Delta_beta   = 0.06;

  // THIS LINE CALCULATES THE NUMBER OF ITERATIONS FOR THE USER'S REFERENCE
  int Ntotal_iterations, counter, cuenta;
  Ntotal_iterations = (int) ( ((beta_end-beta_ini)/Delta_beta)*((gamma_end-gamma_ini)/Delta_gamma)*((adm_end-adm_ini)/Delta_adm)*((as_end-as_ini)/Delta_as)*((Mdm_end-Mdm_ini)/Delta_Mdm) );
  
  printf("  Running %d iterations in parameter space\n", Ntotal_iterations);
  
  counter=0;
  
  //DEFINITION OF THE VECTOR R[i]

  for(i=0;i<K;i++)
    {
      R[i] = (i+1.0);
      if(i==0)
	R[i] = R[i]*0.1;
    }
  
  for(beta = beta_ini; beta <= beta_end; beta = beta + Delta_beta)
    { 
      for(gamma = gamma_ini; gamma <= gamma_end;  gamma = gamma + Delta_gamma)
	{
	  for(a_dm = adm_ini; a_dm <= adm_end;  a_dm = a_dm + Delta_adm)
	    {     
	      for(a_s = as_ini; a_s <= as_end;  a_s = a_s + Delta_as)
		{
		  // THIS IS TO AVOID INDETERMINATIONS IN SOME DENOMINATORS
		  if (a_s != a_dm)
		    {
		      for(M_dm = Mdm_ini; M_dm < Mdm_end;  M_dm = M_dm + Delta_Mdm)
			{
			  for(i = 0; i < K; i ++)
			    {
			      evaluate_integral(R[i], beta, gamma, a_dm, a_s, M_dm, M_s, &sig_p[i], &R_arcmin[i]);
			      //printf("%d R = %lf I=%lf\n",i, log10(R[i]), sig_p[i]);
			    }
			  
			  // ALLOC MEMORY FOR THE INTERPOLATION IN EACH STEP
			  gsl_interp_accel *acc  =  gsl_interp_accel_alloc ();
			  gsl_spline *spline     =  gsl_spline_alloc (gsl_interp_cspline, 90);
			  
			  gsl_spline_init (spline, R_arcmin, sig_p, 90);
			  
			  chi2 = 0.0;
			  for(j=0; j<NUMBER_OF_SIGMA_BINS; j++)
			    {
			      sig_p_i = gsl_spline_eval (spline, R_o[j], acc);
			      chi2 += ( (sig_p_i-sig_p_o[j])*(sig_p_i-sig_p_o[j]) ) / numb_error[j];
			      //printf("%16.8lf\t %16.8lf\n", sig_p_i, chi2);
			    }
			  
			  // FREE MEMORY
			  gsl_spline_free (spline);
			  gsl_interp_accel_free (acc);
			  
			  // CALCULATION OF XI SQUARE TO FIND THE OPTIMIZED PARAMETERS
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
			  
			}// for Mdm	      
		    }// if for the singularity
		}//for as
	    }// for adm
	}//for gamma
    }//for beta
  
  //printf ("\n");
  
  best_model=fopen("best_model.dat","w"); 
  
  for(i = 0; i < K; i ++)          
    {      
      evaluate_integral(R[i], beta_min, gamma_min, a_dm_min, a_s_min, M_dm_min, M_s_min, &sig_p[i], &R_arcmin[i]);
      fprintf(best_model,"%16.8e %16.8e\n", R_arcmin[i], sig_p[i]); 
    }
  
  fclose(best_model);
  
  // THIS ARE THE OPTIMIZED PARAMETERS
  
  printf("xi=%16.8lf beta=%16.8lf a_dm/pc=%16.8lf a_s/pc=%16.8lf M_dm/Msun=%16.8elf M_s/Msun=%16.8e gamma=%lf\n",
	 chi_min, beta_min, a_dm_min, a_s_min, (1.0e5)*M_dm_min, (1.0e5)*M_s_min, gamma_min);
  
  script = fopen( "script.gpl", "w" );
  fprintf( script, "plot 'sigma.dat'pt 7 ps 1.2\n" );
  //fprintf( script, "replot 'best_model.dat' u 1:2 w l\n");
  //fprintf( script, "reset\n");
  fprintf(script, "set grid\nset terminal png\nset output 'sigma.png'\nset nokey\n");
  fprintf( script, "set title 'Sigma Proyectada vs R'\n" );
  fprintf( script, "set xrange [0:40]\n" );
  fprintf( script, "set yrange [4:14]\n" );
  fprintf( script, "set xlabel 'Radius in Arcmin'\n" );
  fprintf( script, "set ylabel 'Projected velocity dispersion in km/s'\n" );
  //fprintf( script, "plot 'sigma.dat'pt 7 ps 1.2\n" );
  fprintf( script, "replot 'best_model.dat' u 1:2 w l\n");
  //fprintf( script, "pause -1\n");
  fclose(script);
  
  warn = system("gnuplot script.gpl");
  
  return(warn);
}
