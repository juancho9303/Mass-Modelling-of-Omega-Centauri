import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#CARGO EL ARCHIVO CON LOS DATOS DE LA DISPERSION DE VELOCIDADES PROYECTADA      
dis,sig = np.loadtxt("sigma.dat",unpack='False')    
for i in range(len(dis)):
  if sig[i] == np.min(sig):
    print sig[i], dis[i], i
    
def func(x,a,b,c):
  return a*x + b*x**2 + c
"""    

#for R in range(0.01, 20, 0.01)

np.piecewise(R/a_s, [R/a_s < 1-1e-9, R/a_s >= 1+e-9], [X_s: (1.0 / sqrt(1.0 - (R/a_s)**2)) * np.log((1.0+(np.sqrt(1.0-(R/a_s)**2)))/(R/a_s)), X_s: (1.0 / np.sqrt((R/a_s)**2 - 1.0)) * acos(1.0/(R/a_s))])
   
def I_R()
  y = (M_s/(2.0*np.pi*a_s**2*GAMMA*(1.0-(R/a_s)**2)*(1.0-(R/a_s)**2)))*((2.0+(R/a_s)**2)*X_s-3.0);
  return y
   
def delta(r,fdfasfsaf)
  y = G/(np.pi*I_R*GAMMA)
  return y
    
def alpha()
  y = ((1.0 - beta*((R*R)/(r*r)))*((r)/(np.sqrt(r*r-R*R))))/2.0*np.pi
  return y
    
def A():
  y = M_s**2*a_s*(12*(a_s + r)**4*np.log10((a_s + r)/(r)) - a_s*(25*a_s**3 + 52*a_s**2*r + 42*a_s*r**2 + 12*r**3))/(12* a_s**5*(a_s + r)**4)
  return y

def D():
  y = M_dm**2*a_dm*(12*(a_dm + r)**4*np.log10((a_dm + r)/(r)) - a_dm*(25*a_dm**3 + 52*a_dm**2*r + 42*a_dm*r**2 + 12*r**3))/(12* a_dm**5*(a_dm + r)**4)
  return y

def B():
  a1 = 2.0*a_s**3*(a_s-a_dm)**4*a_dm*a_dm*(a_s+r)**2*(a_dm+r)
  a2 = 2.0*(a_s-a_dm)**4*(a_s+r)**2*(a_dm+r)*np.log10(r)
  a3 = 2.0*a_dm**2*(6.0*a_s**2-4.0*a_s*a_dm+a_dm*a_dm)*(a_s+r)*(a_s+r)*(a_dm+r)*np.log10(a_s+r)
  a4 = 2.0*a_s**4+4.0*a_s**3*r-2.0*a_dm**2*r*(a_dm+r)+3.0*a_s*a_dm*(-a_dm*a_dm+a_dm*r+2.0*r*r)+a_s*a_s*(7.0*a_dm*a_dm+7.0*a_dm*r+2.0*r*r)
  a5 = 2.0*a_s*a_s*(a_s-4.0*a_dm)*(a_s+r)*(a_s+r)*(a_dm+r)*np.log10(a_dm+r)
  y = M_s*M_dm*((-1.0)/(a1))*(a2-a3+a_s*((a_s-a_dm)*a_dm*a4-a5));
  return y

def C():
  b1 = 2.0*a_dm**3*(a_dm-a_s)**4*a_s*a_s*(a_dm+r)**2*(a_s+r)
  b2 = 2.0*(a_dm-a_s)**4*(a_dm+r)**2*(a_s+r)*np.log10(r)
  b3 = 2.0*a_s**2*(6.0*a_dm**2-4.0*a_dm*a_s+a_s*a_s)*(a_dm+r)*(a_dm+r)*(a_s+r)*np.log10(a_dm+r)
  b4 = 2.0*a_dm**4+4.0*a_dm**3*r-2.0*a_s**2*r*(a_s+r)+3.0*a_dm*a_s*(-a_s*a_s+a_s*r+2.0*r*r)+a_dm*a_dm*(7.0*a_s*a_s+7.0*a_s*r+2.0*r*r)
  b5 = 2.0*a_dm*a_dm*(a_dm-4.0*a_s)*(a_dm+r)*(a_dm+r)*(a_s+r)*np.log10(a_s+r)
  y = M_dm*M_s*((-1.0)/(b1))*(b2-b3+a_dm*((a_dm-a_s)*a_s*b4-b5))
  return y
  
def int(r,fsafdsa)
  y = alpha*(A+B+C+D)
  return y

#PROGRAMA PARA INTEGRAR ------------------------------------------------

def H5(x):
  y=120*x-160*(x**3)+32*(x**5)   #H5
  return y
def func(x):    #integrando
  y=H5(x)**2*exp(-x*x)
  return y
def trap(a,b,n):   #a,b,  limites de integracion
  h=(b-a)/(n-1)
  suma=0.5*(func(a)+func(b))    #primero y ultimos
  for i in arange (1,n-1):
    suma=suma+func(a+i*h)    #da la integral
    return suma*h
  I=trap(-7.0,7.0,1000)
  print(I)
  print(2**5*120.0*sqrt(pi))
#--------------------------------------------------------------------------
  
def sigma()
  y = np.sqrt(delta*int)
  return y


"""

popt, pcov = curve_fit(func, dis, sig)
print np.min(func(dis,popt[0],popt[1],popt[2])), popt[0], popt[1], popt[2]
plt.plot(dis,sig,"o",c = "blue")
for i in range(3):
  plt.plot(dis,func(dis,popt[0],popt[1],popt[2]), "r-", c = "red")
plt.grid()
plt.show()
