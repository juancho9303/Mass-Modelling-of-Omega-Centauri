#include <stdio.h>
#include<stdlib.h>
#include<math.h>

main()
{

 float d[31],temporal;
 FILE *archivo;
 FILE *orden;
 int i,j,p,ID[31],temp;
 int Nd=31;
 int m;

 
 archivo=fopen("distancias.dat","r");
 orden=fopen("ordenados.dat","w");

 for(p=0;p<30;p++)
   {
     fscanf(archivo,"%d %f",&ID[p],&d[p]);
   }

 //sort bubble to distance

 for(i=0;i<(Nd-1);++i)
    {
      for(j=0;j < Nd - 1 - i; ++j)
	{
	 
	    
	      if(d[j] > d[j+1])
		{
		
		
		
		  //	fprintf(orden,"\n");
		temporal = d[j+1];  
		d[j+1]=d[j];          
		d[j]=temporal; 
	      
		temp=ID[j+1];
		ID[j+1]=ID[j];
		ID[j]=temp;   
	    }
	       
	}
      fprintf(orden,"%d %f\n",ID[j],d[j]);
    }

 fclose(archivo);
 fclose(orden);
}
	    


