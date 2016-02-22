#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#define a 5.0
#define pi 3.141592
#define M 5000
#define gamma 1.2
 
  double int_1 (double r, void * params) {
//  double a = *(double *) params;
  double int_1 = 1.0/(r*pow((r+a),5));
  return int_1;
  }

int
main (void)
{  

  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (2000);
  
  double result, error;
  double expected = -4.0;
//  double a = 1.0;

  gsl_function F;
  F.function = &int_1;
//  F.params = &a;

  gsl_integration_qags (&F, 1, 1000, 0, 1e-4, 2000,
                        w, &result, &error); 

  printf ("result          = % .18f\n", result);
  printf ("exact result    = % .18f\n", expected);
  printf ("estimated error = % .18f\n", error);
  printf ("actual error    = % .18f\n", result - expected);
  printf ("intervals       = %zu\n", w->size);

  gsl_integration_workspace_free (w);

  return 0;
}
