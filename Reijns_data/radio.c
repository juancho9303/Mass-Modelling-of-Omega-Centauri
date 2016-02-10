#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 1966

int main (void)
{
  int i;

  double a[N],b[N],c[N],d[N],e[N],f[N],g[N],h[N],j[N],k[N],l[N],AR[N],DEC[N],AR_c[N],DEC_c[N],dis[N];
  
  FILE *data, *center, *radius, *script;
  data = fopen("table4.dat", "r");
  center = fopen("center.dat", "r");
  radius = fopen("distancia.dat" , "w");
  
  for(i=0; i<N; i++)
    {  
      fscanf(data, "%lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\n", &a[i], &b[i], &c[i], &d[i], &e[i], &f[i], &g[i], &h[i], &j[i], &k[i], &l[i], &AR[i], &DEC[i]);
      
      if ( 120.0 < j[i] && j[i] < 400.0 )
	{
	  fscanf(center, "%lf\t %lf\n", &AR_c[i], &DEC_c[i]);
	  dis[i] = sqrt(pow((AR[i]-AR_c[i]),2) + pow((DEC[i]-DEC_c[i]),2));
	  fprintf(radius,"%lf\t %lf\t %lf\n", dis[i], j[i], k[i]);
	}      
    }
    
  fclose(data);
  fclose(center);
  fclose(radius); 
  
  script = fopen( "script.gpl", "w" );
  fprintf(script, "set grid\nset terminal png\nset output 'velocity_vs_rad.png'\nset nokey\n");
  fprintf( script, "set title 'Velocity vs Radius'\n" );
  fprintf( script, "plot 'distancia.dat' u 1:2\n");
  fclose(script);
  
  //  printf("E. Cinetica: %Lf\t E. Potencial: %Lf\t E. Total: %Lf\n", kin_tot, pot_tot, tot_energy);
  
}
