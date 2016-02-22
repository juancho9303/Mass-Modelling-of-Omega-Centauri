#include <stdio.h>
#include <math.h>

#define N 10
#define M 10e-3
#define GAMMA 1.2
#define G 43007.1

  double r, a, a1, a2, a3, a4, a5, b1, b2, b3, b4, b5, beta, a_s, a_dm, M_dm, M_s, M_dmM_s, M_sM_dm, M_dmM_dm, M_sM_s, tot;

int main(){
  
  a = 1.0;
  a_dm = 1.0;
  a_s = 1.1;
  M_s = 2.0;
  M_dm = 2.0;
  beta = 2.0;
  
  FILE *int1, *int2, *int3, *int4, *int5, *graf1, *graf2, *graf3, *graf4, *graf5; 
  
  int1=fopen("integrando1.dat","w");
  int2=fopen("integrando2.dat","w");
  int3=fopen("integrando3.dat","w");
  int4=fopen("integrando4.dat","w");
  //int5=fopen("integrando5.dat","w");
  
  for(r = 0.1; r < 100; r=r+0.01){
    
    M_sM_s = M_s*M_s*a_s*(((12.0*pow((a_s+r),4.0)*log((a_s+r)/(r))) - a_s*(25.0*a_s*a_s*a_s+52.0*a_s*a_s*r+42.0*a_s*r*r+12.0*r*r*r))/(12.0*pow(a_s,5)*pow((a_s+r),4)));
    
    fprintf(int1,"%16.8e\t %16.8e\n", r, M_sM_s);
    
    a1 = 2.0*pow(a_s,3.0)*pow((a_s-a_dm),4.0)*a_dm*a_dm*(a_s+r)*(a_s+r)*(a_dm+r);
    a2 = 2.0*pow((a_s-a_dm),4.0)*(a_s+r)*(a_s+r)*(a_dm+r)*log(r);
    a3 = 2.0*a_dm*a_dm*(6.0*a_s*a_s-4.0*a_s*a_dm+a_dm*a_dm)*(a_s+r)*(a_s+r)*(a_dm+r)*log(a_s+r);
    a4 = 2.0*pow(a_s,4.0)+4.0*pow(a_s,3.0)*r-2.0*a_dm*a_dm*r*(a_dm+r)+3.0*a_s*a_dm*(-a_dm*a_dm+a_dm*r+2.0*r*r)+a_s*a_s*(7.0*a_dm*a_dm+7.0*a_dm*r+2.0*r*r);
    a5 = 2.0*a_s*a_s*(a_s-4.0*a_dm)*(a_s+r)*(a_s+r)*(a_dm+r)*log(a_dm+r);
    
    M_sM_dm = M_s*M_dm*a_dm*((-1.0)/(a1))*(a2-a3+a_s*((a_s-a_dm)*a_dm*a4-a5));
    //printf ("M_sM_dm  = %16.8e\t % .18f\n", r, M_sM_dm);
    fprintf(int2,"%16.8e\t %16.8e\n", r, M_sM_dm);
    
    b1 = 2.0*pow(a_dm,3)*pow((a_dm-a_s),4.0)*a_s*a_s*(a_dm+r)*(a_dm+r)*(a_s+r);
    b2 = 2.0*pow((a_dm-a_s),4.0)*(a_dm+r)*(a_dm+r)*(a_s+r)*log(r);
    b3 = 2.0*a_s*a_s*(6.0*a_dm*a_dm-4.0*a_dm*a_s+a_s*a_s)*(a_dm+r)*(a_dm+r)*(a_s+r)*log(a_dm+r);
    b4 = 2.0*pow(a_dm,4.0)+4.0*pow(a_dm,3.0)*r-2.0*a_s*a_s*r*(a_s+r)+3.0*a_dm*a_s*(-a_s*a_s+a_s*r+2.0*r*r)+a_dm*a_dm*(7.0*a_s*a_s+7.0*a_s*r+2.0*r*r);
    b5 = 2.0*a_dm*a_dm*(a_dm-4.0*a_s)*(a_dm+r)*(a_dm+r)*(a_s+r)*log(a_s+r);
    
    M_dmM_s = M_dm*M_s*a_s*((-1.0)/(b1))*(b2-b3+a_dm*((a_dm-a_s)*a_s*b4-b5));
    //printf ("M_dmM_s  = %16.8e\t % .18f\n", r, M_dmM_s);
    fprintf(int3,"%16.8e\t %16.8e\n", r, M_dmM_s);
    
    M_dmM_dm = M_dm*M_dm*a_dm*(((12.0*pow((a_dm+r),4.0)*log((a_dm+r)/(r))) - a_dm*(25.0*a_dm*a_dm*a_dm+52.0*a_dm*a_dm*r+42.0*a_dm*r*r+12.0*r*r*r))/(12.0*pow(a_dm,5)*pow((a_dm+r),4)));
  
  fprintf(int4,"%16.8e\t %16.8e\n", r, M_dmM_dm);
  
  //tot = M_sM_s + M_dmM_dm + M_sM_dm + M_dmM_s;
  //fprintf(int5,"%16.8e\t %16.8e\n", r, tot);
      
  }
  
  fclose(int1);
  fclose(int2);
  fclose(int3);
  fclose(int4);
  //fclose(int5);
  
/*  
  graf1=fopen("grafica1.gpl","w");
  fprintf(graf1, "set terminal png\nset output 'integrando1.png'\nset nokey\n");
  fprintf( graf1, "set title 'Sigma Proyectada vs R'\n" );
  fprintf( graf1, "plot 'integrando1.dat' u 1:2 w l\n");
  fclose(graf1);
  
  graf2=fopen("grafica2.gpl","w");
  fprintf(graf2, "set terminal png\nset output 'integrando2.png'\nset nokey\n");
  fprintf( graf2, "set title 'Sigma Proyectada vs R'\n" );
  fprintf( graf2, "plot 'integrando2.dat' u 1:2 w l\n");
  fclose(graf2);
  
  graf3=fopen("grafica3.gpl","w");
  fprintf(graf3, "set terminal png\nset output 'integrando3.png'\nset nokey\n");
  fprintf( graf3, "set title 'Sigma Proyectada vs R'\n" );
  fprintf( graf3, "plot 'integrando3.dat' u 1:2 w l\n");
  fclose(graf3);
  
  graf4=fopen("grafica4.gpl","w");
  fprintf(graf4, "set terminal png\nset output 'integrando4.png'\nset nokey\n");
  fprintf( graf4, "set title 'Sigma Proyectada vs R'\n" );
  fprintf( graf4, "plot 'integrando4.dat' u 1:2 w l\n");
  fclose(graf4);
*/
  graf1=fopen("grafica1.gpl","w");
  fprintf( graf1, "plot 'integrando1.dat' u 1:2 w l\n");
  fprintf( graf1, "replot 'integrando2.dat' u 1:2 w l\n");
  fprintf( graf1, "replot 'integrando3.dat' u 1:2 w l\n"); 
  fprintf(graf1, "set terminal png\nset output 'integrando1.png'\nset nokey\n");
  fprintf( graf1, "set title 'Sigma Proyectada vs R'\n" );
  fprintf( graf1, "replot 'integrando4.dat' u 1:2 w l\n");
  
  fclose(graf1);
  
  /*
  graf2=fopen("grafica2.gpl","w");
  fprintf(graf2, "set terminal png\nset output 'integrando2.png'\nset no	key\n");
  fprintf( graf2, "set title 'Sigma Proyectada vs R'\n" );
  fprintf( graf2, "plot 'integrando5.dat' u 1:2 w l\n");
  
  fclose(graf2);
  */
  
  system("gnuplot grafica1.gpl");
  //system("gnuplot grafica2.gpl");
  //system("gnuplot grafica3.gpl");
  //system("gnuplot grafica4.gpl");
    
}
