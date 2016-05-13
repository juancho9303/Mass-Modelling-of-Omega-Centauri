#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#define N 7344

int main (void)
{
  int i, j, warn, counter;
  double lim[15];
  lim[0]=0.01;
  //lim[1]=0.02;
  lim[1]=0.03;
  //lim[3]=0.04;
  lim[2]=0.05;
  //lim[5]=0.09;
  lim[3]=0.11;
  //lim[7]=0.14;
  lim[4]=0.17;
  //lim[9]=0.208;
  lim[5]=0.26;
  //lim[11]=0.32;
  lim[6]=0.4;
  //lim[13]=0.51;
  lim[7]=0.64;
  //lim[15]=0.8;
  lim[8]=1.0;
  //lim[17]=1.25;
  lim[9]=1.56;
  //lim[19]=1.95;
  lim[10]=2.44;
  //lim[21]=3.05;
  lim[11]=3.81;
  //lim[23]=4.76;
  lim[12]=5.96;
  //lim[25]=7.45;
  lim[13]=9.31;
  //lim[27]=11.64;
  lim[14]=14.55;
  //lim[29]=18.18;
  counter = 0.0;
  
  double pos[N],vel[N],err[N],vel_sq, sum;
  double mean, median, skew, sd, kurtosis, error, sum1, sum2;
  
    
  FILE *velo, *sig, *script;
  velo = fopen("organized.dat", "r");
  sig = fopen("sigma.dat","w");
  
  for(i=0.0; i<N; i++)
    {        
      fscanf(velo, "%lf %lf %lf", &pos[i], &vel[i], &err[i]);
    } 
      
  //printf ("The sample mean is %g\n", mean); 
  fclose(velo);
  
  for(j=0.0; j<14.0; j=j+1.0)
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
	      error = 1.0/(sqrt(counter)) * pow(sum1,-0.5) * sum2;

	    }
	}
	
      printf("Particulas entre %lf y %lf es %d\n",lim[j],lim[j+1], counter);
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
      fprintf( script, "set logscale x\n" );
      fprintf( script, "set xlabel 'Radius in Arcmin'\n" );
      fprintf( script, "set xrange [0.01:20]\n" );
      //fprintf( script, "set logscale x 10\n" );
      fprintf( script, "set ylabel 'Sigma (km/s)'\n" );
      fprintf( script, "plot 'sigma.dat' u 1:2:3 pt 7 ps 1 with errorbars\n");
      fclose(script);
      
      warn = system("gnuplot script1.gpl");
      
      return(warn);
      

}



