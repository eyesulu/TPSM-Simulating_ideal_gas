# Project B: Thermodynamics Snookered

## Overview
This project simulates a basic thermodynamic system (ideal gas) using a ball simulation model in a container. The code was created on Spyder 5.1.5 and is structured into the following modules:

## Modules

### 1. `ball` module:
This module contains the `Ball` class, which represents a ball with mass, position, velocity, and radius.

#### Class: Ball
- **Initialization**: To initialize the class, specify the mass, position vector, velocity vector, and radius of the ball.
- **Example parameters**:
  - `mass = 1`
  - `position = [1,1]`
  - `velocity = [0,0]`
  - `radius = 0.5`
  
  Example code to create a ball:
  ```python
  ball = Ball(1, [1, 1], [0, 0], 0.5) 

- **Testing**: All methods in the `Ball` class are demonstrated in the `testing_ball` module.

### 2. `simulation` module:
This module contains functions for simulating the system, including:
- A function for finding minimum values in an array and their indices.
- A function to generate a Maxwell-Boltzmann distribution.
- The `Simulation` class for running the thermodynamic simulation.

#### Class: Simulation
- **Initialization**: To initialize the simulation, specify the number of balls, radius of balls, radius of the container, mass of balls, initial temperature, and whether to display the initial state of the simulation.
- **Example parameters**:
  - `number of balls = 20`
  - `radius of balls = 0.5`
  - `radius of container = 10`
  - `mass of balls = 1`
  - `initial temperature = 2`
  - `initial temperature = 2`
  
  Example code to create a simulation:
  ```python
  sim = Simulation(20, 0.5, 10, 1, 2, show=True)

- **Testing**: All methods in the `Simulation` class are demonstrated in the `testing_simulation` module.

### 3. `testing_ball` module:
This module contains tests for the `ball` module to ensure the methods are working correctly.
- **Note**: Run the `ball` module before running the `testing_ball` module to ensure that the class is available for testing.

### 4. `testing_simulation` module:
This module contains tests for the `simulation` module to ensure its methods are working correctly, and it includes the plotting of graphs.
- **Note**: Run the `ball` and `simulation` modules before running the testing_simulation module.
- **Running Instructions**:
  - Run each cell separately.
  - If a cell takes too long to run, try adjusting the number of frames or iterations.
 
## Dependencies




	








