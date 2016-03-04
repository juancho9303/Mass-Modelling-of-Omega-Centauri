import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
     
dis, bri = np.loadtxt("data_ready.dat", usecols=(0,1), unpack= True)  

def I_R(R,Ie):
  
  y = Ie * np.exp(-7.669*((R/288.0)**0.25-1.0))
  return y

popt, pcov = curve_fit(I_R, dis, bri)

print "-------------------------------------------------------------------"
print "Minimo de la funcion  a_fit   GAMMA_Fit  M_Fit"
print "-------------------------------------------------------------------"
print np.min(I_R(dis,popt[0])), popt[0]
print "-------------------------------------------------------------------"

plt.plot(dis,bri,"o",c = "blue")
for i in range(len(bri)):
  plt.plot(dis,(I_R(dis,popt[0])), "r-", c = "red")
plt.grid()
plt.show()