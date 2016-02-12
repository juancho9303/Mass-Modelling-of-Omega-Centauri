#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 1910

int main (void)
{
  int i, warn, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12;
  
  double dis[N], vel[N], e[N], vel_sq1, sum1, sig1, vel_sq2, sum2, sig2, vel_sq3, sum3, sig3, vel_sq4, sum4, sig4, vel_sq5, sum5, sig5, vel_sq6, sum6, sig6, vel_sq7, sum7, sig7, vel_sq8, sum8, sig8, vel_sq9, sum9, sig9, vel_sq10, sum10, sig10, vel_sq11, sum11, sig11, vel_sq12, sum12, sig12;
  
  FILE *velo, *sigma, *script;
  velo = fopen("velocity.dat", "r");
  //  sigma = fopen("sigma.dat" , "w");
  
  for(i=0; i<N; i++)
    {  
      fscanf(velo, "%lf\t %lf\t %lf\n", &dis[i], &vel[i], &e[i]);      
      
      if ( 0.0 < dis[i] && dis[i] < 1.0 )
	{
	  b1++;
	  vel_sq1 = vel[i]*vel[i];
	  sum1 = vel_sq1 + vel_sq1;
	  sig1 = sqrt(sum1/b1);  
	}
      
      if ( 1.0 < dis[i] && dis[i] < 2.0 )
	{
	  b2++;
	  vel_sq2 = vel[i]*vel[i];
	  sum2 = vel_sq2 + vel_sq2;
	  sig2 = sqrt(sum2/b2);
	}
      
      if ( 2.0 < dis[i] && dis[i] < 3.0 )
	{
	  b3++;
	}
      
      if ( 3.0 < dis[i] && dis[i] < 4.0 )
	{
	  b4++;
	}
      
      if ( 4.0 < dis[i] && dis[i] < 6.0 )
	{
	  b5++;
	}
      
      if ( 6.0 < dis[i] && dis[i] < 8.0 )
	{
	  b6++;
	}
      
      if ( 8.0 < dis[i] && dis[i] < 12.0 )
	{
	  b7++;
	}
      
      if ( 12.0 < dis[i] && dis[i] < 16.0 )
	{
	  b8++;
	}
      
      if ( 16.0 < dis[i] && dis[i] < 20.0 )
	{
	  b9++;
	}
      
      if ( 20.0 < dis[i] && dis[i] < 24.0 )
	{
	  b10++;
	}
      if ( 24.0 < dis[i] && dis[i] < 28.0 )
	{
	  b11++;
	}
      if ( 28.0 < dis[i] && dis[i] < 45.0 )
	{
	  b12++;
	}
    }
  printf(" Bin 1: %d\n Bin 2: %d\n Bin 3: %d\n Bin 4: %d\n Bin 5: %d\n Bin 6: %d\n Bin 7: %d\n Bin 8: %d\n Bin 9: %d\n Bin 10: %d\n Bin 11: %d\n Bin 12: %d\n %lf\n", b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, sig1);
  
  fclose(velo);
  //  fclose(sigma); 
  /*  
      script = fopen( "script1.gpl", "w" );
      fprintf(script, "set grid\nset terminal png\nset output 'Sigma_vs_rad.png'\nset nokey\n");
      fprintf( script, "set title 'Sigma vs Radius'\n" );
      fprintf( script, "set xlabel 'Radius in Arcmin'\n" );
      fprintf( script, "set yrange [25:100]\n" );
      fprintf( script, "set ylabel 'Sigma (km/s)'\n" );
      fprintf( script, "plot 'sigma.dat' u 1:2:3 w yerrorbars lt 12, '' using 1:2 w p pt 7 ps 0.8 lt 12\n");
      fclose(script);
      
      warn = system("gnuplot script1.gpl");
      
      return(warn);
  */  
}
