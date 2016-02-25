#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

// SE UTILIZA UN MODELO DE HERNQUIST MODIFICADO PARA ENCONTRAR LA FORMA DE LA DISPERSION DE VELOCIDADES PROYECTADA Y PODER AJUSTAR LOS PARAMETROS QUE MEJOR SE AJUSTEN A LAS OBSERVACIONES

#define K 0.09           //90 parsecs que esta por encima de los ~60pc de omega centauri
#define GAMMA 1.0        //razon masa luminosidad
#define G 43007.1        //unidades de Gadget
#define EPSILON 1e-5     //para suavizar denominadores para evitar indeterminaciones

/* UNIDADES:
 
  Hubble (internal units) = 0.1
  G (internal units) = 43007.1
  UnitMass_in_g = 1.989e+43 
  UnitTime_in_s = 3.08568e+16 
  UnitVelocity_in_cm_per_s = 100000 
  UnitDensity_in_cgs = 6.76991e-22 
  UnitEnergy_in_cgs = 1.989e+53
  
  Unidad de Masa = 10^10 Msun/h
  Unidad de longitud = 1kpc/h
  Unidad de tiempo ~ 1Gyr
*/

struct param{
  double beta, a_s, a_dm, M_s, M_dm;  
};

double R;
int i,j;

// ACA SE DEFINE LA INTEGRAL QUE CONTIENE TODOS LOS TERMINOS DE LAS CONTRIBUCIONES DE MASA ESTELAR Y MATERIA OSCURA 

double integrando(double r, void * params) 
{
  struct param parameters = *(struct param *) params;
  
  double a_s    =  parameters.a_s;
  double a_dm   =  parameters.a_dm;
  double M_s    =  parameters.M_s;
  double M_dm   =  parameters.M_dm;
  double beta   =  parameters.beta;
  double alpha, a1, a2, a3, a4, a5, b1, b2, b3, b4, b5, M_sM_dm, M_dmM_s, M_sM_s, M_dmM_dm;
  
  alpha = ((1.0 - beta*((R*R)/(r*r)))*((r)/(sqrt(r*r-R*R))));
  
  M_sM_s = M_s*M_s*a_s*(((12.0*pow((a_s+r),4.0)*log((a_s+r)/(r))) - a_s*(25.0*a_s*a_s*a_s+52.0*a_s*a_s*r+42.0*a_s*r*r+12.0*r*r*r))/(12.0*pow(a_s,5)*pow((a_s+r),4)));
  
  M_dmM_dm = M_dm*M_dm*a_dm*(((12.0*pow((a_dm+r),4.0)*log((a_dm+r)/(r))) - a_dm*(25.0*a_dm*a_dm*a_dm+52.0*a_dm*a_dm*r+42.0*a_dm*r*r+12.0*r*r*r))/(12.0*pow(a_dm,5)*pow((a_dm+r),4)));
  
  a1 = 2.0*pow(a_s,3.0)*pow((a_s-a_dm),4.0)*a_dm*a_dm*(a_s+r)*(a_s+r)*(a_dm+r);
  a2 = 2.0*pow((a_s-a_dm),4.0)*(a_s+r)*(a_s+r)*(a_dm+r)*log(r);
  a3 = 2.0*a_dm*a_dm*(6.0*a_s*a_s-4.0*a_s*a_dm+a_dm*a_dm)*(a_s+r)*(a_s+r)*(a_dm+r)*log(a_s+r);
  a4 = 2.0*pow(a_s,4.0)+4.0*pow(a_s,3.0)*r-2.0*a_dm*a_dm*r*(a_dm+r)+3.0*a_s*a_dm*(-a_dm*a_dm+a_dm*r+2.0*r*r)+a_s*a_s*(7.0*a_dm*a_dm+7.0*a_dm*r+2.0*r*r);
  a5 = 2.0*a_s*a_s*(a_s-4.0*a_dm)*(a_s+r)*(a_s+r)*(a_dm+r)*log(a_dm+r);
  
  M_sM_dm = M_s*M_dm*((-1.0)/(a1))*(a2-a3+a_s*((a_s-a_dm)*a_dm*a4-a5));
  
  b1 = 2.0*pow(a_dm,3)*pow((a_dm-a_s),4.0)*a_s*a_s*(a_dm+r)*(a_dm+r)*(a_s+r);
  b2 = 2.0*pow((a_dm-a_s),4.0)*(a_dm+r)*(a_dm+r)*(a_s+r)*log(r);
  b3 = 2.0*a_s*a_s*(6.0*a_dm*a_dm-4.0*a_dm*a_s+a_s*a_s)*(a_dm+r)*(a_dm+r)*(a_s+r)*log(a_dm+r);
  b4 = 2.0*pow(a_dm,4.0)+4.0*pow(a_dm,3.0)*r-2.0*a_s*a_s*r*(a_s+r)+3.0*a_dm*a_s*(-a_s*a_s+a_s*r+2.0*r*r)+a_dm*a_dm*(7.0*a_s*a_s+7.0*a_s*r+2.0*r*r);
  b5 = 2.0*a_dm*a_dm*(a_dm-4.0*a_s)*(a_dm+r)*(a_dm+r)*(a_s+r)*log(a_s+r);
  
  M_dmM_s = M_dm*M_s*((-1.0)/(b1))*(b2-b3+a_dm*((a_dm-a_s)*a_s*b4-b5));
  
  double integrando = alpha*(M_sM_s + M_dmM_dm + M_sM_dm + M_dmM_s);
  
  return integrando;
}

int main(){
  
  FILE *mod,*script,*sigma;
  int warn;
  int Nint = 1000; // Numero de intervalos
  double X_s, I_R, s, R, R_min, sig_p, rho, re1, re2, result1, error1, result2, error2, a_dm, a_s, M_s, M_dm, beta, R_o[12], sig_p_o[12], shi2, shi_t;
  struct param params; // structura de parametros
  gsl_integration_workspace *z = gsl_integration_workspace_alloc(Nint);
  gsl_function F1;
  
  sigma=fopen("sigma.dat","r"); 
  
  for(j=0.0; j<12; j++)
    {        
      fscanf(sigma, "%lf %lf", &R_o[j], &sig_p_o[j]);
    }      

  fclose(sigma);
  
  // PARAMETROS
  
  beta  =  0.000001;
  
  mod=fopen("model.dat","w"); 
  
  for(a_dm = 0.005; a_dm <= 0.06;  a_dm = a_dm + 0.01)
    {
      for(a_s = 0.005; a_s <= 0.06;  a_s = a_s + 0.01)
	{
	  if (a_s != a_dm)
	    {
	      for(M_dm = 0.00001; M_dm < 0.0008;  M_dm = M_dm + 0.0002)
		{
		  for(M_s = 0.0001; M_s < 0.0008;  M_s = M_s + 0.0002)
		    {
		      
		      params.beta  =  beta;
		      params.a_s   =  a_s;
		      params.a_dm  =  a_dm;
		      params.M_s   =  M_s;
		      params.M_dm  =  M_dm;       
		      
		      for(R = 0.001; R < K; R = R + 0.001)          
			{                                     
			  
			  F1.function = &integrando;
			  F1.params = &params;
			  
			  gsl_integration_qags(&F1, R, 0.1, 0, 1e-4, Nint, z, &result1, &error1); 
			  //printf ("result  = % .18f\n", result1);
			  
			  re1 = result1;
			  s = R/a_s;
			  
			  if (s < 1.0 - 1.0e-9)
			    {      	
			      X_s = (1.0 / sqrt(1.0 - s*s)) * log((1.0+(sqrt(1.0-s*s)))/(s));
			    }
			  
			  if(s >= 1.0-1.0e-10 && s <= 1.0+1.0e-10)
			    if(s = 1.0)
			      {
				X_s = 1.0;
			      }
			  
			  if (s > 1.0 + 1.0e-9)
			    {	       	  
			      X_s = (1.0 / sqrt(s*s - 1.0)) * acos(1.0/s);
			    }      
			  
			  I_R = (M_s/(2.0*M_PI*a_s*a_s*GAMMA*(1.0-s*s)*(1.0-s*s)))*((2.0+s*s)*X_s-3.0);
			  //Sigma proyectada
			  sig_p = sqrt(((G)/(GAMMA*(I_R)*M_PI))*re1);
			  //Convierto de kiloparsecs a minutos de arco
			  R_min = R/(4.80839)*3437.75;
			  //printf("%16.8e\t %16.8e\n", R, I_R);
			  fprintf(mod,"%16.8e\t %16.8e\n", R_min, sig_p);
			  
			  for(j=0.0; j<12.0; j=j+1.0)
			    {
			      if( R_min - R_o[j] < 0.2 && R_min - R_o[j] > -0.2 )
				{
				  if(sig_p > 0.0)
				  {
				  shi2 = (sig_p-sig_p_o[j])*(sig_p-sig_p_o[j]);
				  shi_t = shi_t + shi2;
				  }
				  //printf("%lf\n", shi_t);       
				}
			    }
			}		  
		      //gsl_integration_workspace_free(z);
		      fprintf(mod,"#%16.8e\t %16.8e\t %16.8e\t %16.8e\n", a_dm, a_s, M_dm, M_s);
		      printf("%16.8e\t %16.8e\t %16.8e\t %16.8e\n", a_dm, a_s, M_dm, M_s);
		      printf("%lf\n", shi_t); 
		      shi_t = 0.0;
		      fprintf(mod,"\n");
		    } 
		}
	    }
	}
    }  
  
  fclose(mod);
  gsl_integration_workspace_free(z);
  
  script = fopen( "script.gpl", "w" );
  fprintf( script, "plot 'sigma.dat'pt 7 ps 1.2\n" );
  fprintf(script, "set grid\nset terminal png\nset output 'sigma.png'\nset nokey\n");
  fprintf( script, "set title 'Sigma Proyectada vs R'\n" );
  fprintf( script, "set xrange [0:40]\n" );
  fprintf( script, "set yrange [4:14]\n" );
  fprintf( script, "set xlabel 'Radius in Arcmin'\n" );
  fprintf( script, "set ylabel 'Projected velocity dispersion in km/s'\n" );
  fprintf( script, "replot 'model.dat' u 1:2 w l\n");
  fclose(script);
  
  warn = system("gnuplot script.gpl");
  
  return(warn);
}
