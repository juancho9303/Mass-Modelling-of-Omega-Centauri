#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 1910
#define ARRAYSIZE 13

int counter(double pos[], int n);
  
int main (void)
{
  int i, j, s, n[12], warn, lim[13];
  lim[0]=0.0;
  lim[1]=1.0;
  lim[2]=2.0;
  lim[3]=3.0;
  lim[4]=4.0;
  lim[5]=6.0;
  lim[6]=8.0;
  lim[7]=12.0;
  lim[8]=16.0;
  lim[9]=20.0;
  lim[10]=24.0;
  lim[11]=28.0;
  lim[12]=45.0;
  
  double pos[N],vel[N],err[N],blah[N];
  
  FILE *velo, *sigma, *script;
  velo = fopen("velocity.dat", "r");
  sigma = fopen("sigma.dat" , "w");
  
  for(i=0; i<N; i++)
    {        
      fscanf(velo, "%lf\t %lf\t %lf\n", &pos[i], &vel[i], &err[i]);
      
      for(j=0; j<12; j=j++)
	{
	  if ( lim[j] < pos[i] && pos[i] < lim[j+1] )
	    {
	      n[j] = counter(pos[i], ARRAYSIZE);
	      printf("%d\n", n[j]);
	      blah[i] = sqrt(vel[i]*vel[i]*0.5);
	      fprintf(sigma,"%lf\t %lf\n", pos[i], blah[j]);
	    }
	}	   
    }
  
  
  fclose(velo);
  fclose(sigma); 
  
  script = fopen( "script.gpl", "w" );
  fprintf(script, "set grid\nset terminal png\nset output 'sigma_vs_rad.png'\nset nokey\n");
  fprintf( script, "set title 'Sigma vs Radius'\n" );
  fprintf( script, "set xlabel 'Radius in Arcmin'\n" );
  //fprintf( script, "set yrange [100:350]\n" );
  fprintf( script, "set ylabel 'Sigma (km/s)'\n" );
  fprintf( script, "plot 'sigma.dat' u 1:2\n");
  fclose(script);
  
  warn = system("gnuplot script.gpl");
  
  return(warn);
  
}

  int counter(double pos[], int n){
    int flag, j, i, lim[13];

    for(j=0; j<12; j=j++)
      {
	if ( lim[j] < pos[i] && pos[i] < lim[j+1] )
	  {
	    flag=flag+1;
	    printf("%d",flag);
	  }
	return flag;
      }
  }

