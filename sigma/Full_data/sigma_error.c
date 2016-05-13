#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#define N 3749

int main (void)
{
  int i, j, warn, lim[13], counter;
  
  lim[0]=0.0;
  lim[1]=1.0;
  lim[2]=2.0;
  lim[3]=3.0;
  lim[4]=4.0;
  lim[5]=7.0;
  lim[6]=10.0;
  lim[7]=13.0;
  lim[8]=16.0;
  lim[9]=20.0;
  lim[10]=24.0;
  lim[11]=28.0;
  lim[12]=45.0;
  counter = 0.0;
  
  double pos[N],vel[N],err[N],vel_sq, sum;
  double mean, median, skew, sd, kurtosis, error, sum1, sum2;
  
    
  FILE *velo, *sig, *script;
  velo = fopen("velocity.dat", "r");
  sig = fopen("sigma.dat","w");
  
  for(i=0.0; i<N; i++)
    {        
      fscanf(velo, "%lf %lf %lf", &pos[i], &vel[i], &err[i]);
    } 
      
  //printf ("The sample mean is %g\n", mean); 
  fclose(velo);
  
  for(j=0.0; j<12.0; j=j+1.0)
    {
      counter = 0.0;
      for(i = 0.0;(i)<=N; i++)
	{
	  if( lim[j] <= pos[i] && pos[i] < lim[(j)+1] )
	    {
	      counter++;
	      gsl_sort (vel, 1, counter);
	      mean = gsl_stats_mean(vel, 1, counter);
	      median = gsl_stats_median_from_sorted_data (vel, 1, counter);
	      sd = gsl_stats_sd (vel, 1, counter);
	      skew = gsl_stats_skew(vel, 1, counter);
	      kurtosis = gsl_stats_kurtosis (vel, 1, counter); 

	      sum1 = (vel[i]-mean)*(vel[i]-mean);
	      sum1++;
	      sum2 = (vel[i]-mean)*err[i];
	      sum2++;
	      //error = 1.0/(sqrt(counter)) * (pow(,(-0.5))) * () * err[i];
	      error = 1.0/(sqrt(counter)) * pow(sum1,-0.5) * sum2;

	    }
	}
	
      printf("Particulas entre %d y %d es %d\n",lim[j],lim[j +1], counter);
      printf("Error %lf\n", error);
      printf ("The sample mean is %g\n", mean);
      printf ("The median is %g\n", median);
      printf ("Standard deviation is %g\n", sd);
      printf ("Skewness is %g\n", skew);
      printf ("kurtosis is %g\n\n", kurtosis);
      //printf("Sigma in the bin is: %.2f\n\n", sig);
      fprintf(sig, "%lf\t %lf\t %lf\t %d\n", (lim[j]+lim[j+1])*0.5, sd, error, counter);
      
    }
    fclose(sig);
    
    script = fopen( "script1.gpl", "w" );
      fprintf(script, "set grid\nset terminal png\nset output 'Sigma_vs_rad.png'\nset nokey\n");
      fprintf( script, "set title 'Sigma vs Radius'\n" );
      //fprintf( script, "set logscale x\n" );
      fprintf( script, "set xlabel 'Radius in Arcmin'\n" );
      //fprintf( script, "set yrange [25:100]\n" );
      fprintf( script, "set ylabel 'Sigma (km/s)'\n" );
      fprintf( script, "plot 'sigma.dat' u 1:2:3 pt 7 ps 1 with errorbars\n");
      fclose(script);
      
      warn = system("gnuplot script1.gpl");
      
      return(warn);
      

}



