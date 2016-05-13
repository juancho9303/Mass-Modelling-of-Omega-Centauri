#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#define N 7343
#define N_bin 490

int main (void)
{
  int i, j, k, warn, counter;
  counter = 0.0;
  
  double pos[N],vel[N],err[N],vel_sq, sum, rmax, rmin, vel_aux[N_bin], err_aux[N_bin];
  double mean, median, skew, sd, kurtosis, error, sum1, sum2;
  
  FILE *velo, *sig, *script;
  velo = fopen("organized.dat", "r");
  sig = fopen("sigma.dat","w");
  
  for(i=0; i<=N; i++)
    {        
      fscanf(velo, "%lf %lf %lf", &pos[i], &vel[i], &err[i]);
    } 
  
  fclose(velo);
  counter = 0.0;
  
  for(i = 0; i<=N; i++)
    {
      vel_aux[i%N_bin] = vel[i];
      err_aux[i%N_bin] = err[i];
      counter++;
      sum1 = 0.0;
      sum2 = 0.0;
      
      if( counter%N_bin==0)
	{
	  rmin   =   rmax;
	  rmax   = pos[i];
	  
	  gsl_sort (vel_aux, 1, N_bin);
	  mean = gsl_stats_mean(vel_aux, 1, N_bin);
	  median = gsl_stats_median_from_sorted_data (vel_aux, 1, N_bin);
	  sd = gsl_stats_sd (vel_aux, 1, N_bin);
	  skew = gsl_stats_skew(vel_aux, 1, N_bin);
	  kurtosis = gsl_stats_kurtosis (vel_aux, 1, N_bin); 
	  
	  //Aca se calcula el error
	  //for(j=0; j<270; j++)
	  //{
	      sum1 = (vel_aux[i%N_bin]-mean)*(vel_aux[i%N_bin]-mean);
	      sum1++;


	      sum2 = (vel_aux[i%N_bin]-mean)*err_aux[i%N_bin];
	      sum2++;


	      error = 1.0/(sqrt(270)) * pow(sum1,-0.5) * sum2;
	      //}
	  
	  printf("El rango del bin es de %lf hasta %lf\n", rmin, rmax);
	  printf("Particulas entre %lf y %lf es %d\n",rmin,rmax, 300);
	  printf("Error %lf\n", error);
	  printf ("The sample mean is %g\n", mean);
	  printf ("The median is %g\n", median);
	  printf ("Standard deviation is %g\n", sd);
	  printf ("Skewness is %g\n", skew);
	  printf ("kurtosis is %g\n\n", kurtosis);
	  fprintf(sig, "%lf\t %lf\t %lf\t %d\n", (rmin+rmax)*0.5, sd, error, 270);
	}
    }
  fclose(sig);
  
  script = fopen( "script1.gpl", "w" );
  fprintf(script, "set grid\nset terminal png\nset output 'Sigma_vs_rad.png'\nset nokey\n");
  fprintf( script, "set title 'Sigma vs Radius'\n" );
  fprintf( script, "set xlabel 'Radius in Arcmin'\n" );
  fprintf( script, "set xrange [0:12]\n" );
  fprintf( script, "set ylabel 'Sigma (km/s)'\n" );
  fprintf( script, "plot 'sigma.dat' u 1:2:3 pt 7 ps 1 with errorbars\n");
  fclose(script);
  
  warn = system("gnuplot script1.gpl");
  
  return(warn);
  
}



