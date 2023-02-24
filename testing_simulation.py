#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 16:48:14 2021
'''
Simulation is a subclass of Ball class
Has following parameters = (number of balls, radius of balls, radius of container, mass of balls, initial temperature, show the initial state)
'''
@author: aisulu
"""

import ball as bl
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import curve_fit
import simulation 
#%%
'''
define parameters
'''
N = 100 #number of balls
rb = 1 #radius of balls
rc = 70 #radius of container
m = 1 #mass of balls
T = 1 #temperature
n = 100 #number of frames 
n_bins = 20 #number of bins for histograms
#%%
'''
show that the simulation __init__ creates N_balls that do not overlap
in case N_balls is too big for the container given radius of balls and container then an Excpetion is raised 
'''
a = Simulation(70, rb, 20, m, T, show=True) 
'''
returns a collection of balls within the container, highest number of balls is 173 if radius of balls is 1 and radius of container is 20
'''
#%%
'''
show that run simulation works without overlapping or escaping balls
'''
b = Simulation(100, rb, 50, m, T, show = False )
b.run(n,True)
'''
returns an animation of balls moving with randomly generated velocities and positions
'''
#%%
'''
plot evolution of kinetic energy
'''
a = Simulation(20, rb, 20, m, T, show=False) 
a.run(1000, False)
a.E_vs_time()
#%%
'''
plot evolution of total momentum 
'''
a.p_vs_time()
#%%
'''
plot evolution of pressure
'''
a.P_vs_time()
#%%
'''
plot the histogram of distances between balls and the origin
'''
c = Simulation(N, rb, rc, m, T, show = False )
c.run(n, False)
c.hist_d_to_center(n_bins)

#plt.savefig ("Histogram of distance to origin")
#%%
'''
plot the histogam of distances between balls
'''
d = Simulation(N, rb, rc, m, T, show = False )
d.run(n, False)
d.hist_d_between_balls(n_bins)
#plt.savefig ("Histogram of distance between balls")
'''
returns a histogram of distances
i expected a uniform distribution but the histogram looks like a normal distribution with mean at radius of container
'''
#%%
'''
plot histogram of speeds of each ball
'''
e = Simulation(1000, rb, 100, m, T, show = False )
e.run(10, False)
e.hist_speeds(n_bins)
#plt.savefig ("Histogram of speeds")
'''
shows a Maxwell-Bolzmann distribution as more balls and frames are used as expected
'''
#%%
'''
plot a graph of pressure versus temperature
'''
temps = np.linspace(1, 15, 30) #temperature array
iters = 5 #number of iterations
P = 0 #initial pressure
T = 0 #initial temperature
pressures = [] #array of pressures 
temps_final =[] #array of final temperatures as it changes 

for i in range (len(temps)):
    for j in range (iters):
        sim = Simulation(20, rb, 30, m, temps[i], show=False)
        sim.run(500)
        P += sim.get_pressure()
        T = sim.get_T()
    pressures.append(P)
    temps_final.append(T)

#fit functions     
fit,cov = np.polyfit(temps_final,pressures,1,cov=True) 
fit_equation = np.poly1d(fit)

#plotting 
plt.plot(temps_final, fit_equation(temps_final), label='The fit fucntion')
plt.plot(temps_final, pressures, '.', label = 'Data from Simulation')

plt.legend()
plt.xlabel("Temperature")
plt.ylabel("Pressure")
plt.title("Pressure versus temperature (500 frames)")
plt.show()

#%%
'''
plot how pressure and temperature changes when radius changes
'''
radius = list(np.arange(0.1, 2, 0.1)) #array of radius
iters = 5 #number of iterations
P2 = 0 #initial pressure
T2 = 0 #initial temperature
pressures2=[] #array of new pressures
temperatures2=[]#array of new temperatures

for i in range (len(radius)):
    for j in range (iters):
        sim = Simulation(20, radius[i], 30, m, 10, show=False)
        sim.run(500)
        P2 += sim.get_pressure()
        T2 = sim.get_T()
    pressures2.append(P2)
    temperatures2.append(T2)

#Plot of pressure versus radius
plt.subplot(1, 2, 1)
vol = (np.array(radius)) #radius

#fit function
fit2,cov2 = np.polyfit(vol ,pressures2,1,cov=True)
fit_equation2 = np.poly1d(fit2)
#plotting
plt.plot (vol, pressures2, '.', label = 'Data from Simulation')
plt.plot(vol, fit_equation2(vol), label = 'The fit function')

plt.legend()
plt.xlabel("Radius")
plt.ylabel("Pressure")
plt.title("Pressure versus  radius (500 frames)")
plt.show()

#Plot of temperature versus radius
plt.subplot(1, 2, 2)
#plotting 
plt.plot (radius, temperatures2, '.', label = 'Data from Simulation')

plt.legend()
plt.xlabel("Radius")
plt.ylabel("Temperature")
plt.title("Temperature versus radius (500 frames)")
plt.show()
#%%
'''
plot how pressure changes with number of particles to fit to van der Waals
'''
n = np.arange(1, 70, 5) #number of particles
iters = 5 #number of iterations 
P3 = 0 #initial pressure
pressures3=[] #pressure obtained from simulation
temp = 20 #temperature
rad = 20 #radius of the container 

for i in range (len(n)):
    for j in range (iters):
        sim = Simulation(n[i], rb, rad, m, temp , show=False)
        sim.run(500)
        P3 += sim.get_pressure()
    pressures3.append(P3)

#%%
#use van der Waals law to fit the data 
def vand(n, a, b): #a and b are constant to be determined
    p = (n* temp -a*((n **2)/rad**1)+a*b*((n **3)/rad**2))/(rad**1 - n*b)
    return p

#use ideal gas law to fit in the data
def ideal_gas (n, a): #a is a constant 
    p = a*n*temp/rad
    return p

initial_guess = [4, 40] #initial guess for van der Waals equation
po,po_cov=sp.optimize.curve_fit(vand,n,pressures3,initial_guess)

initial_guess2 = 0.5 #initial guess for ideal gas law
po2,po_cov2=sp.optimize.curve_fit(ideal_gas,n,pressures3,initial_guess2)

plt.plot(n, vand(np.array(n), po[0], po[1]), label='van der Waals predictions' )
plt.plot (n, ideal_gas(n, po2[0]), label = 'Ideal gas')
plt.plot (n, pressures3, 'x', label = 'Data from Simulation')

plt.legend()
plt.xlabel("Number of particles")
plt.ylabel("Pressure")
plt.title("Pressure versus number of particles (500 frames)")
plt.show()

print ("Constant a = %.3f, constant b = %.3f" % (po[0], po[1]))


















                    