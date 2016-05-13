#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//ESTE PROGRAMA USA UNA VERSION SIMPLE DE UNA TECNICA LAGRANGIANA PARA RESOLVER PROBLEMAS DE HIDRODINAMICA EN FLUIDOS IDEALES SIN VISCOSIDAD
//Laura Arboleda Hernández

#define M_total 1



float distancia(int x1, int y1, int x2, int y2)
{
  float d;
  d=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
  return d;
}

int main()

{
  int N=500;
  int x[N],y[N],vx[N],vy[N];
  FILE *pf,*dist,*resul;
  float d[N],temporal=0;
  float a,h;
  int ID[N],vecino_id[N][32],temp=0;
  float m=(float)1/N,V[N];
  int p,f,i,j,k,g;
  float ww,densidad[N],presion[N];
  float vecinos[N][32],W[32],Sx,Sy,Ex,Ey,ex,ey,sx,sy;
  float veloxW,veloyW,SUMx,SUMy,U[N],ax[N],ay[N];
 
  //LEER EL ARCHIVO

  pf=fopen("datos.dat","r");
  dist=fopen("distancias.dat","w");
  resul=fopen("resultados.dat","w");

  printf("%f\n",m);

  for(i=0;i<N;i++) 
	{    
	  fscanf(pf,"%d %d %d %d",&x[i],&y[i],&vx[i],&vy[i]);
	}

  //CALCULAR LAS DISTANCIAS
	 
  for(k=0;k<N;k++) //para cada particula (N=10 particulas) 
    {
      for (j=0;j<N;j++) 
	{
	  d[j]=distancia(x[k],y[k],x[j],y[j]);
	  ID[j]=j;
	}  
      
      //ALGORITMO BURBUJA
      
      for(i=0;i<(N-1);++i)  //Organiza las distancias en orden decreciente,moviendo su respectivo ID, 
	                    //el cual es un indicador del vecino respecto al cual se calcula la distancia. 
	
	{
	  for(j=0;j < N -1 - i; ++j)
	    {
	      
	      if(d[j] > d[j+1])
		{
		  
		  temporal = d[j+1];  
		  d[j+1]=d[j];          
		  d[j]=temporal; 
		  
		  temp=ID[j+1];
		  ID[j+1]=ID[j];
		  ID[j]=temp;
		  
		}
	      
	    }
	  
	}
      
      //CALCULO DE h
      a=0.0;
      for(g=0;g<N;g++)
	{
	  a=a+d[g];
	}
      h=0.5*(a/N); //(a/N) : distancia promedio
      
      //SELECCION DE LOS 32 VECINOS MAS CERCANOS

      ww=0.0;
      for(p=1;p<=33;p++)
	{
	  vecinos[k][p]=d[p];
	  vecino_id[k][p]=ID[p];
	  W[p]=(1/M_PI*h*h)*exp(-(vecinos[k][p]*vecinos[k][p])/h*h);
	  ww=ww+W[p]; //sumatoria de los W. Todas las particulas tienen la misma masa m, entonces sale de la sumatoria.
	  
	}
      
      //DENSIDAD

      densidad[k]=ww*m;

      //PRESION
      
      presion[k]=0.5*densidad[k];
      
      Sx=0.0;
      Sy=0.0;
      Ex=0.0;
      Ey=0.0;
      SUMx=0.0;
      SUMy=0.0;
    }


  for(k=0;k<N;k++) //para cada particula (N=500 particulas) 
    {  
      for(j=0;j<32;j++) 
	{	  
	  f=vecino_id[k][j]; //va al vecino ESPECIFICO de una particula dada.
	  
	  //terminos para calcular la aceleracion 
	  ex=(x[k]-x[f])*exp(-(vecinos[k][j]*vecinos[k][j])/h*h);
	  ey=(y[k]-y[f])*exp(-(vecinos[k][j]*vecinos[k][j])/h*h);
	  
	  sx=(presion[f]/(densidad[f]*densidad[f]))*ex;
	  sy=(presion[f]/(densidad[f]*densidad[f]))*ey;
	  
	  Sx=Sx+sx;
	  Sy=Sy+sy;
	  
	  Ex=Ex+ex;
	  Ey=Ey+ey;
	  
	  //terminos para calcular la energia interna
	  veloxW=(-2/M_PI*h*h*h*h)*(vx[k]-vx[f])*ex;
	  veloyW=(-2/M_PI*h*h*h*h)*(vy[k]-vy[f])*ey;
	  
	  SUMx=SUMx+veloxW;
	  SUMy=SUMy+veloyW;
	  
	}
      
      //ACELERACION DEL FLUIDO EN CADA PARTICULA
      
      ax[k]=(m*2/(M_PI*h*h*h*h))*(Ex*(presion[k]/(densidad[k]*densidad[k]))+Sx);
      ay[k]=(m*2/(M_PI*h*h*h*h))*(Ey*(presion[k]/(densidad[k]*densidad[k]))+Sy);
      
     
      //ENERGIA ASOCIADA A CADA PARTICULA DEL FLUIDO
      
      U[k]=(presion[k]/(densidad[k]*densidad[k]))*m*(SUMx+SUMy);
      
      //MAGNITUD DE LA VELOCIDAD DE CADA PARTICULA
      
      V[k]=sqrt(vx[k]*vx[k] + vy[k]*vy[k]);
      
      
      //imprime en el archivo distancias.dat los ID y las distancias organizadas
      
		for (i=0;i<(N);++i)
		{
		fprintf(dist,"%d %f\n",ID[i],d[i]);
		}
		fprintf(dist,"\n");
	       
      //Escribimos las posiciones, la magnitud de la velocidad, presion, densidad y energia interna para cada particula
     
      fprintf(resul,"%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n",x[k],y[k],V[k],ax[k],ay[k],presion[k],densidad[k],U[k]);
	
      
    }
 
  
  fclose(pf);
  fclose(dist);
  fclose(resul);
      
  return 0;
  
}
