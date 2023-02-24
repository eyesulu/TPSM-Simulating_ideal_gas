#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 23:32:10 2021

@author: aisulu
This module checks how ball.py works 
All values below are consistent with theoretical predictions
Ball class has following parameters = (mass, position, velocity, radius)
"""
import ball as bl
import numpy as np
#create two balls 
a = bl.Ball(1, [0,0] , [1,1], 1)
b = bl.Ball(1, [6,6] , [-1,-1], 1)

#get position 
print ("Poitions of the balls before the collision are ", a.pos(), " and ", b.pos())
       
#get velocity
print ("Velocities of the balls are " , a.vel()," and " , b.vel())

#find kinetic energy of the system 
KE_before = 0.5*np.dot(a.vel(), a.vel())+0.5*np.dot(b.vel(), b.vel())
print ("Kinetic energy before the collision = %.3f" %KE_before)

#find time to collision
t_to_collision = a.time_to_collision(b)
print ("Time to collision = %.3f" %t_to_collision)

#move balls to theit colllision time
a.move(t_to_collision)
b.move(t_to_collision)

#position of the collision
print ("Poitions of the balls during the collision are ", a.pos()," and ", b.pos())

#perform collision
a.collide(b)
print ("Velocities after the collision of the balls are ", a.vel()," and " , b.vel())

#check if kinetic energy is conserved
KE_after = 0.5*np.dot(a.vel(), a.vel())+0.5*np.dot(b.vel(), b.vel())
print ("KE before = %.3f" % KE_before)
print ("KE after =%.3f " % KE_after)

