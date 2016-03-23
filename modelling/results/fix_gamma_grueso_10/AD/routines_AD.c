
double integrando(double r, void * params) 
{
  struct param parameters = *(struct param *) params;
  
  double a_s    =  parameters.a_s;
  double a_dm   =  parameters.a_dm;
  double M_s    =  parameters.M_s;
  double M_dm   =  parameters.M_dm;
  double beta   =  parameters.beta;
  double R      =  parameters.R;

  double alpha, b1, b2, b3, b4, b5, c1, c2, c3, c4, c5;
  double M_dmM_s, M_sM_s, M_dmM_dm, M_sM_dm, aux1, aux2;
  
  alpha = (1.0 - (beta*R*R)/(r*r))  * (r/sqrt(r*r-R*R));
  aux1  = a_s+r;
  aux2  = a_dm+r;
  
  // THIS IS A IN OUR NOTATION
  M_sM_s = M_s*M_s*a_s* ( -(25.0*a_s*a_s*a_s + 52.0*a_s*a_s*r + 42.0*a_s*r*r + 12.0*r*r*r)/(12.0*pow(a_s,4)*pow(aux1,4)) +  log(aux1/r)/pow(a_s,5.0));
  
  // THIS IS D IN OUR NOTATION
  M_dmM_dm = M_dm*M_dm*a_dm* ( -(25.0*a_dm*a_dm*a_dm + 52.0*a_dm*a_dm*r + 42.0*a_dm*r*r + 12.0*r*r*r)/(12.0*pow(a_dm,4)*pow(aux2,4)) +  log(aux2/r)/pow(a_dm,5.0));
  
    double integrando = alpha*(M_sM_s + M_dmM_dm + M_sM_dm + M_dmM_s);
  
  return integrando;
}
