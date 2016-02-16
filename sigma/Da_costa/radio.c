#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 1355

int main (void)
{
  int i, warn;

  double dis[N],vel[N];
  
  FILE *data, *velo, *script;
  data = fopen("data.dat", "r");
  velo = fopen("velocity.dat" , "w");
  
  for(i=0; i<N; i++)
    {  
      fscanf(data, "%lf\t %lf\n", &dis[i], &vel[i]);
      
      if ( 120.0 < vel[i] && vel[i] < 400.0 )
	{
	  fprintf(velo,"%lf\t %lf\n", dis[i], vel[i]);
	}      
    }
    
  fclose(data);
  fclose(velo); 
  
  script = fopen( "script.gpl", "w" );
  fprintf(script, "set grid\nset terminal png\nset output 'velocity_vs_rad.png'\nset nokey\n");
  fprintf( script, "set title 'Velocity vs Radius'\n" );
  fprintf( script, "set xlabel 'Radius in Arcmin'\n" );
  //fprintf( script, "set yrange [100:350]\n" );
  fprintf( script, "set ylabel 'Velocity (km/s)'\n" );
  fprintf( script, "plot 'velocity.dat' u 1:2\n");
  fclose(script);
  
  warn = system("gnuplot script.gpl");
  
  return(warn);
  
}
