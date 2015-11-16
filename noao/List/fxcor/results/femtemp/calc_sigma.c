#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 63
#define pi 3.141592
#define cons 0.33
#define G 6.67e-11 
#define fwhma 4.566  //hay que definir el valor real de esta mierda

int main (void)
{  
  FILE *file_fwhm, *sig, *script;
  file_fwhm = fopen("fwhm.txt","r");
  sig = fopen("sigma.txt","w");
  script = fopen( "script_sigma.gpl", "w" );
  int i;
  long double fwhm[N], fwhmi[N], sigma[N], mass[N], r_e[N];
  
  for (i = 1; i < N; i++) 
    { 	
      fscanf(file_fwhm, "%Lf\t %Lf\n", &fwhm[i], &r_e[i]);  
      //printf("%Lf\t %Lf\n", fwhm[i], r_e[i]);
      fwhmi[i] = sqrt((fwhm[i]*fwhm[i])-(fwhma*fwhma));
      sigma[i] = fwhmi[i]/2.35;
      mass[i] = (r_e[i]*sigma[i]*sigma[i])/(cons*G);
      fprintf(sig, "%Lf\t %Lf\n", sigma[i], mass[i]);
      //printf("%Lf\t %Lf\n", sigma[i], mass[i]);
    }  
    
  fclose(sig);  
  fclose(file_fwhm);
  fprintf(script, "set terminal png\nset output 'sigma.png'\nset nokey\n");
  fprintf( script, "set title 'Ga'\n" );
  fprintf( script, "set size square\n");
  fprintf( script, "set logscale y\n");
  fprintf( script, "plot 'sigma.txt' u 2 w d\n");
  fclose(script);
  system("gnuplot script_sigma.gpl");  
  return 0;  
}
