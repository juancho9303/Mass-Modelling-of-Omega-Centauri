#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#define a 5.0
#define pi 3.141592
#define M 5000
#define gamma 1.2

//Variables Globales
//double s;


double X_s (double s, void * params) {
//  double a = *(double *) params;
  double X_s = 1.0 / sqrt(1.0 - s*s) * log(1.0+(sqrt(1-s*s))/(s));
  return X_s;
}
/*
  double I_R () {
  double I_R = 2.0 + X_s;
  return I_R;
  }
*/  
  double I_s (double s, void * params) {
//  double a = *(double *) params;
  double I_s = M/(2*pi*a*a*gamma*(1-s*s)*(1-s*s))*((2+s*s)*(1.0 / sqrt(1.0 - s*s) * log(1.0+(sqrt(1-s*s))/(s)))-3);
  return I_s;
  }
//struct X_s
//{
//  double s;
//};

int
main (void)
{  

  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (100);
  
  double result, error;
  double expected = -4.0;
//  double a = 1.0;

  gsl_function F;
  F.function = &I_s;
//  F.params = &a;

  gsl_integration_qags (&F, 0, 1, 0, 1e-5, 100,
                        w, &result, &error); 

  printf ("result          = % .18f\n", result);
  printf ("exact result    = % .18f\n", expected);
  printf ("estimated error = % .18f\n", error);
  printf ("actual error    = % .18f\n", result - expected);
  printf ("intervals       = %zu\n", w->size);

  gsl_integration_workspace_free (w);

  return 0;
}