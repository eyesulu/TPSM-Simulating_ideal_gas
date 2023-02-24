#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 16:40:05 2021

@author: aisulu
"""

import numpy as np
import pylab as pl 
from itertools import combinations
import random
import warnings
import matplotlib.pyplot as plt
import statistics


'''
Ball class

Initialised as Ball (mass, posi, velocity, rad) where 
    mass = mass of the ball/particle
    posi = position of the ball/particle, an array with x and y components
    velocity = velocity of the ball/particle, an array with x and y components
    rad = radius of the ball 

Has following methods:
    pos() - returns the position of the ball
    vel() - returns the velocity of the ball
    mass() - returns the mass of the ball
    move(dt) - moves the ball by delta time
    time_to_collision(other ball) - returns time to collision with 
            the other ball
        the time the method returns is the minimum positive time 
        if this is not possible (it is possible to have complex results or 
                                    only negative times)
        then it returns 1e6 to avoid problems in Simulation 
    collide(other ball) - performs the collision and updates the velocities
        of each ball
    patch() - returns patch of the ball
'''
class Ball:
    
    def __init__(self, mass, posi, velocity, rad):
        
        #setting attributes
        self._m = mass
        self._posi = np.array(posi)
        self._vel = np.array(velocity)
        self._rad = rad
        
        #setting patches, works with the container too
        if self._rad == 1:
            self._patch = pl.Circle(self._posi, self._rad , fc='r')
        else:
            self._patch=pl.Circle(self._posi, self._rad, edgecolor = 'black',\
                                  fill = False, ls='solid' ) 
            
    def pos(self):
        return self._posi
    
    def vel(self):
        return self._vel
    
    def move(self, dt):
        self._posi = self._posi+self._vel * dt
        self._patch.center = self._posi #updates the center of the patch 
        return self
    
    def mass(self):
        return self._m
    
    def time_to_collision (self, other):
        
        r = self._posi - other._posi #position vector between two objects 
        V = self._vel - other._vel #difference in velocities
        
        #equation is different for ball-ball and ball-container collision 
        #therefore need to account for that
        if self._rad == other._rad:
            R = self._rad + other._rad
        else: 
            R = self._rad - other._rad
        
        #if two objects are stationary returns 1e6
        if np.dot(V, V) == 0:
            t1=1e6
            return t1
        
        #divide the equation from script by components
        first = 2*np.dot(r, V) #2r.V
        sec = first**2 #(2r.V)^2
        thi = np.dot(r, r) -R**2 #r^2 - R^2
        four = np.dot(V, V) #V^2
        
        #set up conditions for which value of time to return
        if (sec-4*four*thi)<0: #if complex roots
            return 1e6
        #rounding solutions to 4 decimal places to avoid problems 
        t1 = np.around (((-first+(sec-4*four*thi)**0.5)/(2*four)), 4)
        t2 = np.around (((-first-(sec-4*four*thi)**0.5)/(2*four)), 4)
        
        if (sec-4*four*thi)==0: #if discriminant is equal to zero
            if t1<0: #if the only solution is negative
                return 1e6
            else: #if the solution is 0 or positive
                return t1
        elif (sec-4*four*thi)>0: #if discriminant is positive
            if t1<=0: #if first solution is negative
                if t2<=0: #if second solution is also negative
                    return 1e6
                else: #if the second solution is zero or positive 
                    return t2 
            elif t2<=0:# same logic as before but with the second solution
                if t1<=0:
                    return 1e6
                else:
                    return t1
            else: #if both solutions are positive 
                return (min(abs(t1), abs(t2)))
            
    def collide (self,other):
        
        rnotnorm = (self._posi - other._posi ) #not normalised position vector
        r = rnotnorm / (np.sqrt (rnotnorm[0]**2+ rnotnorm[1]**2)) # normalising
        v1 = np.dot (self._vel, r) #projection ie parallel component to r
        v2 = np.dot (other._vel, r) # parallel component to r of the second ball
        x = np.array([1,0]) #x unit vector
        y = np.array ([0,1]) #y unit vector
        
        if self._rad ==other._rad: #if collisions between two balls
            v1new = v2 * r#equation for v1'
            v2new = v1 * r  # equation for v2'
            v1p = -v1*r +self._vel#perpendicular component which is unchanged
            v2p = -v2*r +other._vel#perpendicular component of the second ball
            
            #express our parallel and perpendicular components 
            #in terms of x and y
            x1 = np.dot(v1new, x) + np.dot(v1p, x) 
            y1 = np.dot (v1new, y) + np.dot(v1p, y)
            x2 = np.dot(v2new,x) +  np.dot(v2p, x)
            y2 = np.dot (v2new, y) + np.dot(v2p, y)
            
            #updating velocities and rounding to 5 decimal places
            self._vel = np.around(np.array ([x1, y1]), 5)
            other._vel =np.around( np.array ([x2, y2]), 5)
            return self, other 
        else:  #if collision with the container using same logic
            v1new = -v1*r
            v1p = -v1*r +self._vel#perpendicular component which is unchanged
            x1 = np.dot (v1new, x) + np.dot(v1p, x)
            y1 = np.dot (v1new, y) + np.dot(v1p, y)
            x2 = 0
            y2 = 0
            self._vel = np.around (np.array ([x1, y1]), 5)
            other._vel = np.around (np.array ([x2, y2]), 5)
            return self, other
        
    def get_patch(self):
        return self._patch



        
            
        
        
    
        