import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
     
dis, bri = np.loadtxt("datos2_ready_parsec_cut.dat", usecols=(0,1), unpack= True)  

#for i in range(len(dis)):
  #if bri[i] == np.min(bri):
    #print bri[i], dis[i], i
    
    #A=100000
def I_R(R,a,A):
    y = np.piecewise(R, 
		   [R < a, R == a, R > a], 
		   [lambda R:  ( ((A*a*a)/((a*a - R*R)*(a*a - R*R))) * (( (2.0 + (R/a)**2) * np.log((1.0 + np.sqrt(1.0 - (R/a)**2))/(R/a)) / np.sqrt(1.0 - (R/a)**2) ) - 3.0) ),
                   lambda R: ( (2.0*A)/(7.5*a*a) ),
                   lambda R: ( ((A*a*a)/((a*a - R*R)*(a*a - R*R))) * (( (2.0 + (R/a)**2) * np.arccos(a/R) / np.sqrt((R/a)**2 - 1.0) ) - 3.0) ) ]) 

popt, pcov = curve_fit(I_R, dis, (bri))
#popt, pcov = curve_fit(I_R, dis, np.log(bri))

print "-------------------------------------------------------------------"
print "Minimo de la funcion  a_fit   GAMMA_Fit  M_Fit"
print "-------------------------------------------------------------------"
print np.min( I_R(dis,popt[0],popt[1]) ), popt[0], popt[1]
print "-------------------------------------------------------------------"

plt.plot(dis, (bri),"o",c = "blue", label = "Data")
plt.plot(dis, (I_R(dis,popt[0],popt[1])) , "r-", c = "red", label = "Fit")    
plt.legend()  
plt.ylabel("I(r)")
plt.xlabel("r (Parcsec)")
plt.grid()
plt.show()



"""
def I_R(R,a,A):
  
  if R < a : 
    y = np.log( ((A*a*a)/((a*a - R*R)*(a*a - R*R))) * (( (2.0 + (R/a)**2) * np.log((1.0 + np.sqrt(1.0 - (R/a)**2))/(R/a)) / np.sqrt(1.0 - (R/a)**2) ) - 3.0) )
  if R == a :
    y = np.log( (2.0*A)/(7.5*a*a) )
  if R > a :
    y = np.log( ((A*a*a)/((a*a - R*R)*(a*a - R*R))) * (( (2.0 + (R/a)**2) * np.arccos(a/R) / np.sqrt((R/a)**2 - 1.0) ) - 3.0) )

  return y

"""








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