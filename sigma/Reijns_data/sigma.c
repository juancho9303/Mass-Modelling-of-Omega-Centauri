#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int function (int, int, double[]);
//int prom();

int main(){
  
  int i, q, w, t, j=0, n[1593];
  double sigma[10], v[10];
  
  for (i=0; i<=8; i=i+2){
    n[j] = function (i, i+2, v);
    j=j+1;
  }
  for (q=8; q<=28; q=q+4){
    n[j] = function (q, q+4, v);
    j=j+1;
  }
  for (w=28; w<=44; w=w+16){
    n[j] = function (w, w+16, v);
    j=j+1;
  }
  
  for (t=0; t<j; t=t+1){
    sigma[t] = sqrt((pow(v[t], 2))/n[t]);
  }
  return (0);
  system("PAUSE");
}

int function (int m, int n, double v[]){
  int i;
  double e, aux;
  FILE *velo, *sigma;
  velo = fopen("velocity.dat", "r");
  sigma = fopen("sigma.dat" , "w");
  int flag;
  double x;
  for (i=0; i<1593; i++) v[i] = 0;
  while (!feof (velo)){
    fscanf(velo, "%lf\t %lf\t %lf\n", &x, &aux, &e);
    if (x >= m && x <= n){
      v[i] += aux;
      flag=flag+1; 
      i=i+1;
    }
    v[i] += v[i]/flag;
  }
  return flag;
}