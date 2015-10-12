//Juan Manuel Espejo Salcedo
//Proyecto sobre condiciones iniciales en Galaxias Esfericas
//Profesor Juan Carlos Muñoz
//Fisica y Astrofisica Computacional

/* -----------------------------------------------------------------------

Este proyecto consiste en generar una estructura esferica que cumpla con unas condiciones iniciales que le permitan cierta estabilidad conforme transcurre el tiempo, es decir, se generan las coordenadas y velocidades de las particulas en el sistema para que este sea un sistema esfericamente estable y que si se llegase a correr en el tiempo, muestre esa estabilidad.

 Se sigue el modelo de Plummer que sugiere ciertas condiciones de densidad y un potencial adecuado para la estructura. 

Se verifica que el sistema es estable ya que cumple con el teorema del virial y se hace un desplazamiento del centro de masa para comodidad. Sin embargo, esto no afecta las posiciones relativas entre las particulas.
 ------------------------------------------------------------------------ */

//Incluyo las librerias de gsl y para las operaciones metematicas

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

//Defino las constantes que voy a usar para las operaciones y ademas un factor de escala que sale del manejo adecuado de las magnitudes

#define N 20001
#define pi 3.141592
#define sca 16.0 / 3*pi

//Inicializo la rutina principal

int main (void)
{
  
  //Genero los archivos donde voy a escribir las coordenadas y velocidades de las particulas que siguen el modelo de Plummer, y defino mis variables de tipo double para tener una alta precision en los calculos
  
  FILE *part, *script;
  part=fopen("particulas.dat","w");
  int i,j;
  double cumulative_mass, ra, posx, posy, posz, xx, yy, vel, velx, vely, velz, pos_cm_x[N], pos_cm_y[N], pos_cm_z[N],cm, velo_x[N], velo_y[N], velo_z[N], p_cm_x, p_cm_y,p_cm_z,v_cm_x,v_cm_y,v_cm_z, mass, theta, phi;
  
  //Llamo librerias y rutinas de GSL para generar numeros aleatorios pues necesito crearlos para las posiciones y velocidades,  
  
  const gsl_rng_type * T;
  gsl_rng * r;
  
  //Utilizo la semilla default para crear los numeros aleatorios
  
  gsl_rng_env_setup();  
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  //Creo las particulas en forma aleaoria en coordenadas esfericas, la coordenada phi va a ser simplemente un numero aleatorio entre 0 y 2*pi, mientras que la otra coordenada angular se halla con el coseno inverso de un numero aleatorio entre -1 y 1
  
  for (i = 1; i < N; i++) 
    {
      mass = 1.0 / N; 
      phi = gsl_rng_uniform (r) * (2*pi);
      theta = acos((gsl_rng_uniform (r)-0.5)*2);

      //Para la distancia radial no voy a utilizar el radio sino la masa acumulada que sera el parametro de entrada para encontrar las distancias de las particulas, segun las deducciones matematicas del modelo de Plummer este valor aleatorio debe cunplir la siguiente relacion:
      
      do
	cumulative_mass = ((gsl_rng_uniform (r)*(1.0-0.00005)+0.00005)/((mass)+i*(1.0/N)));

      //Donde se filtran los resultados que dan negativos y generan indeterminaciones

      while(cumulative_mass > 1.0);

      //la coordenada radial entonces queda en funcion de la masa acumulativa que crece a una unidad de masa a la vez
      
      ra = (1.0 / sqrt(pow(cumulative_mass,(-2.0/3.0)) - 1.0)) / sca; 

      //Paso mis coordenadas a coordenadas cartesianas:
      
      posx = ra * sin(theta) * cos(phi);
      posy = ra * sin(theta) * sin(phi);
      posz = ra * cos(theta);
      
      //Ahora la velocidad, tambien la genero con numeros aleatorios, pero esta vez se debe cumplir una relacion diferente, primero creo numeros aleatorios para xx y yy: 
      
      xx = gsl_rng_uniform (r);
      yy = gsl_rng_uniform (r);

      //Si se cumple la siguiente condicion, entonces genero los numeros aleatorios en cierto rango, para xx entre 0 y 1 y para yy entre 0 y 0.1
      
      if( yy > xx*xx * pow((1.0-xx*xx),3.5) )
	{
	  xx = gsl_rng_uniform (r);
	  yy = gsl_rng_uniform (r)/10.0;
	}
      
      //Genero el valor de la velocidad de acuerdo a la condicion anterior y de nuevo utilizo coordenadas cartesianas: 

      vel = xx * sqrt(2.0) * pow(( 1.0 + ra*ra),-0.25)/pow(sca,(1/3));
      velx = vel * sin( theta ) * cos( phi );
      vely = vel * sin( theta ) * sin( phi );
      velz = vel * cos( theta );

      //Ahora, voy a calcular el centro de masa de la distribucion de particulas, para hacer el debido corrimiento tanto para la posicion como para la velocidad. la idea es que el centro de masa de la distribucion se desplace junto con la estructura esferica
      
      p_cm_x = 0;
      p_cm_y = 0;
      p_cm_z = 0;
      v_cm_x = 0;
      v_cm_y = 0;
      v_cm_z = 0;

	  pos_cm_x[i] = mass * posx;
	  p_cm_x = p_cm_x + pos_cm_x[i];
	  pos_cm_y[i] = mass * posy;
	  p_cm_y = p_cm_y + pos_cm_y[i];
	  pos_cm_z[i] = mass * posz;
	  p_cm_z = p_cm_z + pos_cm_z[i];
	  velo_x[i] = mass * velx;
	  v_cm_x = v_cm_x + velo_x[i];
	  velo_y[i] = mass * vely;
	  v_cm_y = v_cm_y + velo_y[i];
	  velo_z[i] = mass * velz;
	  v_cm_z = v_cm_z + velo_z[i];

	  posx = posx - p_cm_x;
	  posy = posy - p_cm_y;
	  posz = posz - p_cm_z;
            
	  //printf("%lf\t %lf\t %lf\t %lf\t %lf\t %lf\n", p_cm_x, p_cm_y, p_cm_z, v_cm_x, v_cm_y, v_cm_z);

	  //Imprimo en un archivo la masa, posicion y velocidad en todas las coordenadas para todas las particulas de la distribucion

      fprintf (part, "%.7lf\t %.5lf\t %.5lf\t %.5lf\t %.5lf\t %.5lf\t %.5lf\t %.5lf\n", mass, posx, posy, posz, ra, velx, vely, velz);
    }  
    
  //Finalmente y con el fin de visualizar la distribucion esferica que cree voy a hacer un script en bash desde mi programa en C donde especifico los parametros para graficar.

  fclose(part);
  gsl_rng_free (r);
  script = fopen( "script.gpl", "w" );
  fprintf(script, "set terminal png\nset output 'particulas.png'\nset nokey\n");
  fprintf( script, "set title 'Galaxia Esferica'\n" );
  fprintf( script, "set size square\n");
  fprintf( script, "splot [-0.33:0.33][-0.33:0.33][-0.18:0.18]'particulas.dat' u 2:3:4 w d\n");
  fclose(script);
  system("gnuplot script.gpl");
  
  return 0;
  
}

//Fin del programa, ahora para hacer el chequeo de que el sistema es correcto segun el modelo de Plummer, voy a hacer chequeos de energia donde espero ver que se cumpla el teorema del virial.
