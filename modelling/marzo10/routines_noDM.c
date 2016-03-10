// ACA SE DEFINE LA INTEGRAL QUE CONTIENE TODOS LOS TERMINOS DE LAS CONTRIBUCIONES DE MASA ESTELAR Y MATERIA OSCURA Y QUE DEPENDE DE r 


double integrando(double r, void * params) 
{
  struct param parameters = *(struct param *) params;
  
  double a    =  parameters.a;
  double M   =  parameters.M;
  double beta   =  parameters.beta;
  double R      =  parameters.R;

  double alpha, b1, b2, b3, b4, b5, c1, c2, c3, c4, c5, M_sM_dm;
  double M_dmM_s, M_sM_s, M_dmM_dm, aux1, aux2;
  
  alpha = (1.0 - (beta*R*R)/(r*r))  * (r/sqrt(r*r-R*R));
  
  double integrando = alpha * (((12.0*pow((a+r),4.0)*log((a+r)/(r))) - a*(25.0*a*a*a+52.0*a*a*r+42.0*a*r*r+12.0*r*r*r)))/(12.0*pow(a,5)*pow((a+r),4));
  
  return integrando;
}
