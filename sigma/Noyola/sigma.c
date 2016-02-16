#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 43
#define AR_c 201.6969999	
#define DEC_c -47.47947222

int main (void)
{
  int i, warn;

  double dis[N],vel[N],min[N],err[N];
  
  FILE *data, *velo, *script;
  data = fopen("data.dat", "r");
  velo = fopen("velocity.dat" , "w");
  
  for(i=0; i<N; i++)
    {  
      fscanf(data, "%lf\t %lf\n", &dis[i], &vel[i]);
      
      //if ( 120.0 < j[i] && j[i] < 400.0 )
	//{
	  min[i] = dis[i]*(1.0/60.0);
	  err[i] = vel[i]*0.05;
	  fprintf(velo,"%lf\t %lf\t %lf\n", min[i], vel[i], err[i]);
	//}      
    }
    
  fclose(data);
  fclose(velo); 
  
  script = fopen( "script.gpl", "w" );
  fprintf(script, "set grid\nset terminal png\nset output 'velocity_vs_rad.png'\nset nokey\n");
  fprintf(script, "set logscale x\n");
  fprintf( script, "set title 'Sigma vs Radius'\n" );
  fprintf( script, "set xlabel 'Radius in Arcmin'\n" );
  //fprintf( script, "set yrange [100:350]\n" );
  fprintf( script, "set ylabel 'Sigma (km/s)'\n" );
  fprintf( script, "plot 'velocity.dat' u 1:2:3 w yerrorbars lt 12, '' using 1:2 w p pt 7 ps 0.8 lt 12\n");
  
  //plot 'velocity.dat' u 1:2 pt 7 ps 1\n");
  fclose(script);
  
  warn = system("gnuplot script.gpl");
  
  return(warn);
  
}
