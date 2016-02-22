import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import arange
from scipy import integrate
from sys import exit

### Parameters
epsilon = 1e-10
G       = 43007.1
GAMMA   = 1.1
beta    = 1.1
#p0=(1.5, 1.7, 1.5, 1.7)
a_s     = 1.7
a_dm    = 1.5
M_s     = 1.7
M_dm    = 1.5
R       = 0.01
#r     = 5.0

dis,sig = np.loadtxt("sigma.dat",unpack='False')    
#for i in range(len(dis)):
#  if sig[i] == np.min(sig):
#    print sig[i], dis[i], i

def I(R,M_s,a_s):
  s  = R/a_s
  if (s>0.0)*(s<1.0-1.0e-9):
    X_s = (1.0 / np.sqrt(1.0 - (s)**2)) * np.log((1.0 + (np.sqrt(1.0-(s)**2)))/(s))
  
  elif s > 1.0 + 1.0e-9:
    X_s = (1.0 /np.sqrt((s)**2 - 1.0)) * np.arccos(1.0/(s))
  
  return (M_s/(2.0*np.pi*a_s**2*GAMMA*(1.0-(s)**2)*(1.0-(s)**2)))*((2.0+(s)**2)*X_s-3.0)


def A(r, M_s, M_dm, a_s, a_dm):
  return M_s**2*a_s*(12*(a_s + r)**4*np.log10((a_s + r)/(r)) - a_s*(25*a_s**3 + 52*a_s**2*r + 42*a_s*r**2 + 12*r**3))/(12* a_s**5*(a_s + r)**4)


def B(r, M_s, M_dm, a_s, a_dm):
  a1 = 2.0*a_s**3*(a_s-a_dm)**4*a_dm*a_dm*(a_s+r)**2*(a_dm+r)
  a2 = 2.0*(a_s-a_dm)**4*(a_s+r)**2*(a_dm+r)*np.log10(r)
  a3 = 2.0*a_dm**2*(6.0*a_s**2-4.0*a_s*a_dm+a_dm*a_dm)*(a_s+r)*(a_s+r)*(a_dm+r)*np.log10(a_s+r)
  a4 = 2.0*a_s**4+4.0*a_s**3*r-2.0*a_dm**2*r*(a_dm+r)+3.0*a_s*a_dm*(-a_dm*a_dm+a_dm*r+2.0*r*r)+a_s*a_s*(7.0*a_dm*a_dm+7.0*a_dm*r+2.0*r*r)
  a5 = 2.0*a_s*a_s*(a_s-4.0*a_dm)*(a_s+r)*(a_s+r)*(a_dm+r)*np.log10(a_dm+r)
  return M_s*M_dm*((-1.0)/(a1))*(a2-a3+a_s*((a_s-a_dm)*a_dm*a4-a5))
    
def C(r, M_s, M_dm, a_s, a_dm): 
  b1 = 2.0*a_dm**3*(a_dm-a_s)**4*a_s*a_s*(a_dm+r)**2*(a_s+r)
  b2 = 2.0*(a_dm-a_s)**4*(a_dm+r)**2*(a_s+r)*np.log10(r)
  b3 = 2.0*a_s**2*(6.0*a_dm**2-4.0*a_dm*a_s+a_s*a_s)*(a_dm+r)*(a_dm+r)*(a_s+r)*np.log10(a_dm+r)
  b4 = 2.0*a_dm**4+4.0*a_dm**3*r-2.0*a_s**2*r*(a_s+r)+3.0*a_dm*a_s*(-a_s*a_s+a_s*r+2.0*r*r)+a_dm*a_dm*(7.0*a_s*a_s+7.0*a_s*r+2.0*r*r)
  b5 = 2.0*a_dm*a_dm*(a_dm-4.0*a_s)*(a_dm+r)*(a_dm+r)*(a_s+r)*np.log10(a_s+r)
  return M_dm*M_s*((-1.0)/(b1))*(b2-b3+a_dm*((a_dm-a_s)*a_s*b4-b5))
  
def D(r, M_s, M_dm, a_s, a_dm):
  return M_dm**2*a_dm*(12*(a_dm + r)**4*np.log10((a_dm + r)/(r)) - a_dm*(25*a_dm**3 + 52*a_dm**2*r + 42*a_dm*r**2 + 12*r**3))/(12* a_dm**5*(a_dm + r)**4)
    
def argumento(r, R, M_s, M_dm, a_s, a_dm):
  alpha = ((1.0 - beta*((R*R)/(r*r)))*((r)/(np.sqrt(r*r-R*R+epsilon))))  
  return (A(r, M_s, M_dm, a_s, a_dm) + B(r, M_s, M_dm, a_s, a_dm) + C(r, M_s, M_dm, a_s, a_dm) + D(r, M_s, M_dm, a_s, a_dm))*alpha

for r in np.arange(0.1,10,0.1):
  print r,argumento(r, R, M_s, M_dm, a_s, a_dm)

def integral(R, M_s, M_dm, a_s, a_dm):
  return integrate.quad(argumento, R, np.inf, args=(R, M_s, M_dm, a_s, a_dm))

for R in np.arange(0.1,1,0.1):
  print integral(R, M_s, M_dm, a_s, a_dm)

#def integral(R, M_s, M_dm, a_s, a_dm):
#  return integrate.quad(argumento, R, np.inf, args=(R, M_s, M_dm, a_s, a_dm))
  #print integral(R, M_s, M_dm, a_s, a_dm)
  
"""  f1=open('./integral.dat', 'w+')  
  f1.write(' %s. %s\n' % (R, integral(R, M_s, M_dm, a_s, a_dm)))
  f1.write('\n')
  
# print integral(R, M_s, M_dm, a_s, a_dm)

##  R = 30
##  print argumento(r, R, M_s, M_dm, a_s, a_dm)


integral = np.vectorize(integral)
  
def sigma2_p(R, M_s, M_dm, a_s, a_dm):
  return integral(R, M_s, M_dm, a_s, a_dm)*G/(np.pi*I(R,M_s,a_s)*GAMMA) 

popt, pcov = curve_fit(sigma2_p, dis, sig, p0=p0)
print np.min(sigma2_p(dis,popt[0],popt[1],popt[2], popt[3])), popt[0], popt[1], popt[2],popt[3]
plt.plot(dis,sig,"o",c = "blue")
for i in range(3):
  plt.plot(dis,sigma2_p(dis,popt[0],popt[1],popt[2],popt[3]), "r-", c = "red")
plt.grid()
plt.show()


#CARGO EL ARCHIVO CON LOS DATOS DE LA DISPERSION DE VELOCIDADES PROYECTADA  

dis,sig = np.loadtxt("sigma.dat",unpack='False')    
for i in range(len(dis)):
  if sig[i] == np.min(sig):
    print dis[i], sig[i], i
    
r     = np.linspace(0.01,20,0.01)

for R in np.arange(0.01,20,0.01):
  s     = R/a_s
  #X_s = np.piecewise(s, [s < 1.0-1e-9, s >= 1.0+1e-9], [lambda x: (1.0 / sqrt(1.0 - (x)**2)) * np.log((1.0+(np.sqrt(1.0-(x)**2)))/(x)), lambda x: (1.0 /np.sqrt((x)**2 - 1.0)) * np.arccos(1.0/(x))])
  
  #def X_s(s):
      #X_s = np.zeros(np.shape(s))
  if s < 1.0 - 1.0e-9:
    X_s = (1.0 / np.sqrt(1.0 - (s)**2)) * np.log((1.0 + (np.sqrt(1.0-(s)**2)))/(s))
    
  if s > 1.0 + 1.0e-9:
    X_s = (1.0 /np.sqrt((s)**2 - 1.0)) * np.arccos(1.0/(s))
      
  I_R = (M_s/(2.0*np.pi*a_s**2*GAMMA*(1.0-(s)**2)*(1.0-(s)**2)))*((2.0+(s)**2)*X_s-3.0)
  delta = G/(np.pi*I_R*GAMMA)
def alpha(r):
  return ((1.0 - beta*((R*R)/(r*r)))*((r)/(np.sqrt(r*r-R*R))))/2.0*np.pi
  
  def sigma(r,M_s,M_dm,a_s,a_dm):
    y = np.sqrt(delta*int)
    return y
  popt, pcov = curve_fit(sigma, dis, sig)
  print popt[0], popt[1], popt[2], popt[3]
  
  M_s  = popt[0]
  M_dm = popt[1]
  a_s  = popt[2]
  a_dm = popt[3]
  
#INTEGRANDO  
def integral(r,M_s,M_dm,a_s,a_dm):
  y = alpha*(A + B + C + D)
  return y

#PROGRAMA PARA INTEGRAR ------------------------------------------------

def trap(a,b,n):   #a,b,  limites de integracion
  h=(b-a)/(n-1)
  suma=0.5*(int(a)+int(b))    #primero y ultimos
  for i in arange (1,n-1):
    suma=suma+int(a+i*h)    #da la integral
    return suma*h
  I=trap(R,40.0,1000)
  print(I)

#--------------------------------------------------------------------------
  
def sigma(r,M_s,M_dm,a_s,a_dm):
  y = np.sqrt(delta*int)
  return y

popt, pcov = curve_fit(sigma, dis, sig)
print np.min(sigma(dis,popt[0],popt[1],popt[2])), popt[0], popt[1], popt[2]
plt.plot(dis,sig,"o",c = "blue")
for i in range(3):
  plt.plot(dis,sigma(dis,popt[0],popt[1],popt[2]), "r-", c = "red")
plt.grid()
plt.show()
"""
