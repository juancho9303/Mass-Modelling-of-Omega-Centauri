#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 7611
#define AR_c 265.1754
#define DEC_c -53.67433

int main (void)
{
  int i, warn;

  double vel[N],err[N],AR[N],DEC[N],dis[N];
  
  FILE *data, *velo, *script;
  data = fopen("data.dat", "r");
  velo = fopen("velocity.dat" , "w");
  
  for(i=0; i<N; i++)
    {  
      fscanf(data, "%lf\t %lf\t %lf\t %lf\n", &AR[i], &DEC[i], &vel[i], &err[i]);
      
      if ( -12.0 < vel[i] && vel[i] < 48.0 )
	{
	  dis[i] = (sqrt(pow((AR[i]-AR_c),2) + pow((DEC[i]-DEC_c),2)))*60.0;
	  fprintf(velo,"%lf\t %lf\t %lf\n", dis[i], vel[i], err[i]);
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
  fprintf( script, "plot 'velocity.dat' u 1:2:3 w yerrorbars lt 12, '' using 1:2 w p pt 7 ps 0.8 lt 12\n");
  fclose(script);
  
  warn = system("gnuplot script.gpl");
  
  return(warn);
  
}
