import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
     
dis,bri = np.loadtxt("datos2_ready.dat",unpack='False')    
for i in range(len(dis)):
  if bri[i] == np.min(bri):
    print bri[i], dis[i], i

def I_R(R,a,GAMMA,M):
  
  y = np.piecewise(R, [(R < a), (R = a), (R > a)], [lambda R: M/(2.0*np.pi*GAMMA)*((2.0 + (R/a)**2.0)*(np.log((1.0 + np.sqrt(1.0 - (R/a)**2.0))/(R/a))/(np.sqrt(1.0 - (R/a)**2.0))) - 3.0)/(a**2.0*(1.0 - (R/a)**2.0)**2.0),
						    lambda R: 2.0*M/(7.5*a*a*2.0*np.pi*GAMMA),
						    lambda R: M/(2.0*np.pi*GAMMA)*((2.0 + (R/a)**2.0)*( np.arccos(1.0/(R/a))/ np.sqrt((R/a)**2.0 - 1.0) ) - 3.0)/(a**2.0*(1.0 - ((R/a)**2.0)**2.0)])
  return y

popt, pcov = curve_fit(I_R, dis, bri)
print np.min(func(dis,popt[0],popt[1],popt[2])), popt[0], popt[1], popt[2]
plt.plot(dis,bri,"o",c = "blue")
for i in range(3):
  plt.plot(dis,I_R(dis,popt[0],popt[1],popt[2]), "r-", c = "red")
plt.grid()
plt.show()











"""
def I_R(R,a,GAMMA,M):
    I_R = np.zeros(np.shape(R))
    for i,R in enumerate(R):
      if R<a:
        I_R[i]= M/(2.0*np.pi*GAMMA)*((2.0 + (R/a)**2.0)*(np.log((1.0 + np.sqrt(1.0 - (R/a)**2.0))/(R/a))/(np.sqrt(1.0 - (R/a)**2.0))) - 3.0)/(a**2.0*(1.0 - (R/a)**2.0)**2.0)
        elif R=a:
	  I_R[i]= 2.0*M/(7.5*a*a*2.0*np.pi*GAMMA)
	  else:
	    I_R[i]= M/(2.0*np.pi*GAMMA)*((2.0 + (R/a)**2.0)*( np.arccos(1.0/(R/a))/ np.sqrt((R/a)**2.0 - 1.0) ) - 3.0)/(a**2.0*(1.0 - ((R/a)**2.0)**2.0)
	return(I_R)
"""