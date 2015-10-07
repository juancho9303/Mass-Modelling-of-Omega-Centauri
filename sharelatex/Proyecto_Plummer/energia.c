//Juan Manuel Espejo Salcedo
//Proyecto sobre condiciones iniciales en Galaxias Esfericas
//Profesor Juan Carlos Muñoz
//Fisica y Astrofisica Computacional

/* -----------------------------------------------------------------------

En esta parte se calcula la energia cinetica, potencial y total del sistema para verificar que el sistema cumple con el teorema del virial.
 ------------------------------------------------------------------------ */

//Incluyo las librerias necesarias

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Defino el numero de particulas que usare para hacer el calculo de la energia, notese que este programa es computacionalmente mas pesado que el anterior ya que hay que hacer un calculo para cada particula con respecto a las demas

#define N 3000

/* Para un valor de 3000 el valor que se encuentra para las energias es de 

E. Cinetica: 0.048208	 E. Potencial: -0.102709	 E. Total: -0.054501

*/

//Para hacer los calculos de Energia:

//Inicializo la rutina principal:

int main (void)
{
  
  //Creo las varibles y los archivos de lectura y escritura

  int i,j,k,q;

  //Defino una nueva variable n que es el numero de distancias que se van a calcular en el arreglo, obviamente depende de el numero de particulas que voy a usar para calcular las energias

  int n = ((N*N)-N);
  long double m[N], x[N], y[N], z[N], ra[N], vx[N], vy[N], vz[N], pot_energy[n], tot_energy, vel_square[N], kin_energy[N], kin_tot, pot_tot, ma[n],dis[n];
  
  FILE *datos, *pot, *pote;
  pot = fopen("dist_pot.dat", "w");
  datos = fopen("particulas.dat" , "r");

  //le doy un valor inicial de 0 a la energia cinetica

  kin_tot = 0;
  
  //hago un ciclo en el que leo el archivo de datos y calculo los cuadrados de las velocidades, que es necesario para la enrgia cinetica, ademas, defino la enegia cinetcia total como la suma de todas las energias cineticas de la galaxia

  for(i=0; i<N; i++)
    {  
      fscanf(datos, "%Lf\t %Lf\t %Lf\t %Lf\t %Lf\t %Lf\t %Lf\t %Lf\n", &m[i], &x[i], &y[i], &z[i], &ra[i], &vx[i], &vy[i], &vz[i]);
      vel_square[i] = vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
      kin_energy[i] = 0.5 * m[i] * vel_square[i];
      //printf("%Le\t kin:%Le\n", vel_square[i], kin_energy[i]);

      //Encuentro el valor total de energia cinetica

      kin_tot = kin_tot + kin_energy[i];

      //Este nuevo ciclo es para calcular la energia potencial total del sistema, para esto debo calcular todas las distancias del arreglo, es decir todas las distancias asociadas a cada una de las particulas de la galaxia teniendo en cuenta que no me interesa la distancia entre la particula i y la particula i pues es cero. Esto es facil de quitar por medio de un condicional if.
      
      for(j=0; j<N; j++)
	{
	  fscanf(datos, "%Lf\t %Lf\t %Lf\t %Lf\t %Lf\t %Lf\t %Lf\t %Lf\n", &m[j], &x[j], &y[j], &z[j], &ra[j], &vx[j], &vy[j], &vz[j]);
	  if( i != j)
	  {	      
	    dis[j] = sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]));
	      fprintf(pot,"%Lf\t %Lf\n", m[i], dis[j]);
	  }	
	}      
    }

  //Cierro los archivos de donde lei los datos y donde guarde las distancias para la energia potencial

  fclose(datos);
  fclose(pot);
      
  pote = fopen("dist_pot.dat","r");

  //En este nuevo ciclo escaneo el archivo con las masas y distancias
  
  for(q=0; q<n; q++)
    {
      fscanf(pote, "%Lf %Lf", &ma[q], &dis[q]);
      //printf("%d %Lf %Lf\n", q, ma[q], dis[q]);
    }
  
  pot_tot = 0;

  //Y en este calculo ahora si las energias potenciales y las sumo para encontrar la energia potencial total del sistema de particulas

  for(k=0; k<n; k++)
    {
      //printf("%d %Lf %Lf\n", k, ma[k], dis[k]);
      pot_energy[k] = ma[k]*ma[k] / sqrt(dis[k]*dis[k]);
      pot_tot = pot_tot - 0.5 *  pot_energy[k];	      
    }
  //pot_tot = N*pot_tot;
  tot_energy = kin_tot + pot_tot;  

  //Imprimo los valores de las energias apra corroborar que mi sistema es estable y comple con el modelo de Plummer para galaxias esfericas
  
  printf("E. Cinetica: %Lf\t E. Potencial: %Lf\t E. Total: %Lf\n", kin_tot, pot_tot, tot_energy);
  
  fclose(pote);
  
}

//Fin del programa
