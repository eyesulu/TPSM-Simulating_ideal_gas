#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  5 05:26:19 2021

@author: aisulu
"""

import numpy as np
import pylab as pl 
from itertools import combinations
import random
import warnings
import matplotlib.pyplot as plt
import statistics
from ball import Ball

'''
Simulation class

Initialised as Simulation(N_balls, rad_balls,rad_cont, mass, T, show=False) 
where
    N_balls = number of balls 
    rad_balls = radius of balls
    rad_cont = radius of the container
    mass = mass of the balls
    T = starting temperature
    show = determines whether to show the simulation
    
Has following methods:
    KE() - retuns kinetic energy of the system 
    get_min_time() - returns minimum times for the next collision 
        and its indicies in the list 
    next_collision() - performs the next collision 
        finds min time for the next collision
        moves all balls by that time
        performs the collision (or several)
        it also sums up all momentums and energies that hit the container 
    run(n_frames, animate = False) - runs the animation for 'n_frames' frames
        also returns values of distances between balls and 
            between balls and the origin and speeds of balls
    get_pressure () - returns the pressure on the wall of the container
    hist_d_to_center()-returns a histogram of distances between balls and 
                the origin
    hist_d_between_balls - returns a histogram of distances between balls
    get_speeds() - returns all speeds
    hist_speeds() - returns a histogram of speed
    get_T() - returns temperature
    
'''
# a function to determine several minimas and their indicies in the list
def min_val(a):
    indicies=[]
    minim = min(a)
    for index, element in enumerate (a):
        if minim==element :
            indicies.append(index)
    return minim, indicies

#Maxwell-Boltzmann distribution
def maxwell (speed, T, m): #Maxwell-Boltzman distribution
    fc = 4*np.pi*(speed**2)*((m/(2*np.pi*T))**(3/2))*np.exp((-m*speed**2)\
                                                            /(2*T))
    return fc

class Simulation(Ball):
            
    
    def __init__ (self, N_balls, rad_balls,rad_cont, mass, T, show=False):
        self._N_balls = N_balls #number of balls
        self._rad_cont = rad_cont #radius of container
        #container object
        self._container = Ball(1000, [0,0], [0,0], self._rad_cont) 
        self._rad = rad_balls #radius of balls
        self._mass = mass #mass of balls
        self._time =0 #time since beginning 
        #pressure added due to one ball colliding with the container:
        self._pressure = 0  
        self._T = T #initial temperature

        #creating random veocities for balls
        mean_vel = np.sqrt((2*self._T)/self._mass) #mean speed
        #normal distribution:
        velocities = mean_vel*np.random.normal(size=(self._N_balls, 2)) 
        
        #creating random positions for balls
        allx = np.arange(-(self._rad_cont-0.5), +self._rad_cont-0.5,\
                         2*self._rad+0.5)
        #distance between balls is 0.5
        ally = np.arange((-self._rad_cont-0.5), +self._rad_cont-0.5, \
                         2*self._rad+0.5)
        allpos=[] 
        for i in range(len(allx)):
            for j in range(len(ally)):
                pos = (allx[i], ally[j])
                allpos.append(pos)
        #distance from centres of mass to centre of container:
        distance=[]
        for i in range(len(allpos)):
            d = np.sqrt(np.dot(allpos[i], allpos[i]))
            distance.append(d)
        #all possible positons 
        possible_pos = []
        for i in range (len(distance)):
            if distance[i] < (self._rad_cont-0.5-self._rad): 
                #deleting positions outside of the container
                possible_pos.append(allpos[i])
                continue
        random.shuffle(possible_pos) #randomising positions 
        
        #creating list of balls
        positions=[]
        if N_balls>len(possible_pos): #if number of balls is too high
            raise Exception ('Number of balls is to high given the radius \
                             of the container and individual balls')
        for i in range (N_balls):
            positions.append(possible_pos[i])
        list_balls=[]
        for i in range (len(positions)):
            obj = Ball (self._mass,positions[i], velocities[i], self._rad )
            list_balls.append(obj)
        self._list_balls = list_balls
        
        #showing the initial state of the container and balls
        if show:
            ax = pl.axes(xlim=(-self._rad_cont, self._rad_cont), \
                         ylim=(-self._rad_cont, self._rad_cont))  
            patch1=self._container.get_patch()
            ax.add_artist(patch1) 
            for obj in self._list_balls:
                self._ball = obj
                patch2=self._ball.get_patch()
                ax.add_patch(patch2)
            pl.show()
    
    #returns total kinetic energy of the system
    def KE(self):
            totalE=0
            for obj in self._list_balls:
                E = 0.5*self._mass*np.dot(obj.vel(),obj.vel())
                totalE+=E
            return totalE
    
    #returns minimum time to collision with its indicies
    def get_min_time (self):
        #combinations of pairs of balls
        pairs=list(combinations( self._list_balls, 2)) 
        time1=[] #time to collision between balls
        for ball1, ball2 in pairs:
            time_balls = ball1.time_to_collision(ball2)
            time1.append(time_balls)
        time2=[] #time to collision between ball and the container 
        for obj in self._list_balls:
            time_con = obj.time_to_collision(self._container)
            time2.append(time_con)
        time = list(time1)+list(time2) #list of all times to collision
        min_time, index = min_val(time) #finding minimum
        return min_time, index
    
    #returns total momentum of the system 
    def get_mom(self):
        total_mom =0 
        for obj in self._list_balls:
            mo= self._mass*obj.vel()
            mom = np.dot(mo, mo)**0.5
            total_mom+=mom
        return total_mom
    
    #moves all balls by minimum time, performs the collision and 
    #sums up all momentums exerted on the container
    def next_collision (self):
        min_time, index = self.get_min_time()
        #creating pairs of ball-ball + ball-container
        pairs=list(combinations( self._list_balls, 2)) 
        for obj in self._list_balls:
            paircont = (obj, self._container)
            pairs.append(paircont)
        #movinf objects by min time to collision
        for obj in self._list_balls:
            obj.move(min_time)
            
        pairs_list = list(pairs)
        moms=0 #total momentum on the container
        
        for i in range(len(index)):
            #if collided with the container
            if self._container == pairs_list[index[i]][1]:
                mom = self._mass*(np.dot((pairs_list[index[i]][0].vel()), \
                                         (pairs_list[index[i]][0].vel()))**0.5)
                moms+=mom #total momentum
            #performing the collision
            pairs_list[index[i]][0].collide(pairs_list[index[i]][1]) 
            
        self._moms = moms
        self._time+=min_time #time since the start of the simulation
        
    def run(self, num_frames,animate = False):
        self._num_frames = num_frames #number of frames
        if animate: #adding patches
            pl.figure()
            ax =pl.axes(xlim=(-self._rad_cont, self._rad_cont), \
                        ylim=(-self._rad_cont, self._rad_cont))
            ax.add_artist(self._container.get_patch()) 
            for obj in self._list_balls:
                ax.add_patch(obj.get_patch())
                
        distance_to_center =[]
        distance_between_balls=[]
        #creatng pairs of balls
        pairs = list(combinations( self._list_balls, 2)) 
        speeds=[]
        pressures =[]
        T_dist =[] #distribtion of temperature
        times=[] #array of times
        energies =[] #array of kinetic energies
        mom_dist=[] #distribution of momentums
        
        for i in range(num_frames):
            
            self.next_collision()
            #creating a distribution of temperatures:
            T_dist.append(self.KE()/self._N_balls) 
            #array of pressures:
            pressures.append (self._moms/(self._time*self._rad_cont)) 
            times.append(self._time) #array of times
            energies.append(self.KE()) #array of kinetic energies
            mom_dist.append(self.get_mom()) #array of total momentums
            
            for obj in self._list_balls:  
                speed =np.sqrt(np.dot( obj.vel(),  obj.vel()))
                speeds.append(speed) #array of speeds
            for obj in self._list_balls:
                r = obj.pos()
                d = np.sqrt(np.dot(r, r))
                #array of distances between balls and the origin:
                distance_to_center.append(d) 
            for pair in pairs:
                r1=pair[0].pos()
                r2=pair[1].pos()
                r = r1-r2
                d = np.sqrt(np.dot(r, r))
                #array of distances between balls:
                distance_between_balls.append(d) 
            if animate:
                pl.pause(0.1)
                pl.show()

        self._distance_to_center = distance_to_center
        self._distance_between_balls=distance_between_balls
        self._speeds = speeds 
        self._pressures = pressures
        self._T_dist = T_dist
        self._times = times
        self._energies = energies
        self._mom_dist = mom_dist
        
    def get_pressure(self): #returns final pressure
        self._pressure =  np.mean(self._pressures)
        return (self._pressure)
    
    def get_T(self): #returns final temperature
        return (self._T_dist[-1])
    
    def get_speeds (self): #returns speeds
        return self._speeds
    
    '''
    followig methods return plots
    '''
    #histogram of distances to origin
    def hist_d_to_center(self, n_bins): 
        plt.hist(self._distance_to_center, n_bins)
        plt.title("Distance from centers of mass of balls to origin")
        plt.xlabel("Distance")
        plt.ylabel("Frequency")
        plt.show()
    
    #histogram of distacnes between balls
    def hist_d_between_balls(self, n_bins):
        plt.hist(self._distance_between_balls, n_bins)
        plt.title("Distance from centres of mass of balls (pairs of balls)")
        plt.xlabel("Distance")
        plt.ylabel("Frequency")
        plt.show()
    
    #histogram of speeds 
    def hist_speeds(self, n_bins):
        vals=[]
        for i in range (len(self._speeds)):
            val = 2200*maxwell(self._speeds[i], self._T, self._mass)
            vals.append(val)
        
        plt.hist(self._speeds, n_bins)
        plt.plot(self._speeds, vals, '.', \
                 label ='Maxwell-Boltzman distribution prediction')
        plt.xlabel("Distance")
        plt.ylabel("Frequency")
        plt.legend()
        plt.show()
    
    #evolution of pressure with time
    def P_vs_time(self):
        x = self._times
        y = self._pressures
        plt.plot(x, y, 'x')
        plt.xlabel("Time")
        plt.ylabel("Pressure")
        plt.title ("Pressure vs Time")
        plt.show()
    
    #evolution of kinetic energy with time
    def E_vs_time(self):
        x = self._times
        y = self._energies
        plt.plot(x, y, 'x')
        plt.xlabel("Time")
        plt.ylabel("Kinetic energy of the system")
        plt.title ("Kinetic energy vs Time")
        plt.show()
    
    #evolution of total momentums with time
    def p_vs_time(self):
        x = self._times
        y = self._mom_dist
        plt.plot(x, y, 'x')
        plt.xlabel("Time")
        plt.ylabel("Magnitude of momentum of the system")
        plt.title ("Momentum vs Time")
        plt.show()
