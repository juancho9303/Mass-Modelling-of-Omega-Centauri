#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 947
#define AR_c 201.6969999	
#define DEC_c -47.47947222

int main (void)
{
  int i, warn;

  double AR[N],DEC[N],vel[N],err[N],dis[N];
  
  FILE *data, *velo, *script;
  data = fopen("omega.dat", "r");
  velo = fopen("velocity.dat" , "w");
  
  for(i=0; i<N; i++)
    {  
      fscanf(data, "%lf\t %lf\t %lf\t %lf\n",&AR[i],&DEC[i],&vel[i],&err[i]);
      
      if ( 20.0 < vel[i] && vel[i] < 320.0 )
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
