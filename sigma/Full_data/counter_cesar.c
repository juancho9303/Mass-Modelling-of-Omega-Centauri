#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 1910

int main (void)
{
  double i, j, warn, lim[13], counter;
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
  counter = 0.0;
  
  double pos[N],vel[N],err[N];
  
  FILE *velo;
  velo = fopen("velocity.dat", "r");
  
  for(i=0.0; i<N; i++)
    {        
      fscanf(velo, "%lf %lf %lf", &pos[(int)i], &vel[(int)i], &err[(int)i]);
    }      
  
  fclose(velo);
  
  for(j=0.0; j<12.0; j=j+1.0)
    {
      counter = 0.0;
      for(i = 0.0;((int)i)<=N; i++)
	{
	  if( lim[(int)j] < pos[(int)i] && pos[(int)i] < lim[((int)j)+1] )
	    {
	      counter++;
	    }
	}
      printf("Numero de particulas con posicion entre %.2f y %.2f es %d\n",lim[(int)j],lim[(int)j +1], (int)counter);
    }
}


