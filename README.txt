
CID:01874240
The code was created on Spyder 5.1.5
Project B: Thermodynamics Snookered

Consists of following files:

	1.ball module: -Contains class Ball-To initialise the class specify mass, position vector, velocity vector, and radius. 
-Use mass = 1. 
-Example parameters: mass = 1, position = [1,1], velocity = [0,0], radius = 0.5
-Example code: ball = Ball(1, [1,1], [0,0], 0.5)
-All methods are demonstrated in the testing_ball script.

	2.simulation module:
-Contains a function for finding minimum values in an array and their indicies, a function for Maxwell-Boltzmann distribution, and class Simulation
-To initialise the class specify number of balls, radius of balls, radius of container, mass of balls, initial temperature, and whether you would like to see the initial state of the simulation.
-Example parameters: number of balls = 20, radius of balls = 0.5, radius of container = 10, mass of balls =1, initial temperature = 2, show the initial state of the simulation
-Example code: sim = Simulation(20, 0.5, 10, 1, 2, show = True)
-All methods are shown in the testing_simulation script. 

	3.testing_ball module:
-Contains tests for ball module to ensure its methods are working correctly.
-Run ball module prior to running the testing_ball module. 

	4.testing_simulation module:
-Contains tests for simulation module to ensure its methods are working correctly and plots graphs. 
-Run ball and simulation modules prior to running the testing_simulation module.
-Run each cell separately. 
-If a cell is taking too long to run, change the number of frames and/or number of iterations. 


	






