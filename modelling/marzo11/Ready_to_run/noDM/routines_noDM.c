
double integrando(double r, void * params) 
{
  struct param parameters = *(struct param *) params;
  
  double a      =  parameters.a;
  double M      =  parameters.M;
  double beta   =  parameters.beta;
  double R      =  parameters.R;

  double alpha;

  alpha = (1.0 - (beta*R*R)/(r*r))  * (r/sqrt(r*r-R*R));

  double integrando = alpha * M*M*a* ( -(25.0*a*a*a + 52.0*a*a*r + 42.0*a*r*r + 12.0*r*r*r)/(12.0*pow(a,4)*pow((a+r),4)) +  log((a+r)/r)/pow(a,5.0));
  
  return integrando;

}
