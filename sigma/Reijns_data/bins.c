#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 1593

int funcion (m, n);

int main (void)
{
  int i, warn;
  
  double dis[N], v[N], e[N], sig;
  
  FILE *velo, *sigma, *script;
  velo = fopen("velocity.dat", "r");
  sigma = fopen("sigma.dat" , "w");
  
  int flag=0;
  
  for(i=0; i<N; i++)
    {  
      fscanf(velo, "%lf\t %lf\t %lf\n", &dis[i], &v[i], &e[i]);
      
      if ( 0.0 < dis[i] && dis[i] < 2.0 )
	{
	  flag=flag+1;
	  return flag;
	  prom = (v[i]*v[i])/flag;
	  sig = sqrt(prom);	
	  fprintf(sigma,"%lf\t %lf\t %lf\n", dis[i], sig, e[i]);
	  //printf("%lf\t", sig);
	}
      
      if ( 2.0 < dis[i] && dis[i] < 4.0 )
	{
	  flag++;
	  return flag;
	  sig = sqrt((v[i]*v[i])/flag);
	  fprintf(sigma,"%lf\t %lf\t %lf\n", dis[i], sig, e[i]);
	}
      
      if ( 4.0 < dis[i] && dis[i] < 6.0 )
	{
	  flag++;
	  return flag;
	  sig = sqrt((v[i]*v[i])/flag);
	  fprintf(sigma,"%lf\t %lf\t %lf\n", dis[i], sig, e[i]);
	}
      
      if ( 6.0 < dis[i] && dis[i] < 8.0 )
	{
	  flag++;
	  return flag;
	  sig = sqrt((v[i]*v[i])/flag);
	  fprintf(sigma,"%lf\t %lf\t %lf\n", dis[i], sig, e[i]);
	}
      
      if ( 8.0 < dis[i] && dis[i] < 12.0 )
	{
	  flag++;
	  return flag;
	  sig = sqrt((v[i]*v[i])/flag);
	  fprintf(sigma,"%lf\t %lf\t %lf\n", dis[i], sig, e[i]);
	}
      
      if ( 12.0 < dis[i] && dis[i] < 16.0 )
	{
	  flag++;
	  return flag;
	  sig = sqrt((v[i]*v[i])/flag);
	  fprintf(sigma,"%lf\t %lf\t %lf\n", dis[i], sig, e[i]);
	}
      
      if ( 16.0 < dis[i] && dis[i] < 20.0 )
	{
	  flag++;
	  return flag;
	  sig = sqrt((v[i]*v[i])/flag);
	  fprintf(sigma,"%lf\t %lf\t %lf\n", dis[i], sig, e[i]);
	}
      
      if ( 20.0 < dis[i] && dis[i] < 24.0 )
	{
	  flag++;
	  return flag;
	  sig = sqrt((v[i]*v[i])/flag);
	  fprintf(sigma,"%lf\t %lf\t %lf\n", dis[i], sig, e[i]);
	}
      
      if ( 24.0 < dis[i] && dis[i] < 28.0 )
	{
	  flag++;
	  return flag;
	  sig = sqrt((v[i]*v[i])/flag);
	  fprintf(sigma,"%lf\t %lf\t %lf\n", dis[i], sig, e[i]);
	}
      
      if ( 28.0 < dis[i] && dis[i] < 44.0 )
	{
	  flag++;
	  return flag;
	  sig = sqrt((v[i]*v[i])/flag);
	  fprintf(sigma,"%lf\t %lf\t %lf\n", dis[i], sig, e[i]);
	}
    }
  
  fclose(velo);
  fclose(sigma); 
  
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
  
}

int funcion (m, n);
{
  int flag;
}
