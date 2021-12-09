# Orbital Modeling

### CMSE 202 003 Semester Project
#### Authors:

Ethan Hill, Emma Oswald, Patrick Tutt, Steven Vancamp


# Background:
**Science 2:**

Model and visualize phenomena in physics (such as particle collisions, doppler shift, life cycle of a star, asteroid trajectories, star spectral analysis, gravitational redshift, and many other possibilities).
We chose orbital modeling for our project because we felt that we had a good enough background understanding of the concepts at play and that it would be simple enough to understand yet complicated enough to be interesting.

**Project goal:**

Create a python object that can model and create visualizations for a solar system.
We wrote our classes including methods and 3d graphs to represent a modeled solar system and compare that to real world data to see how well our model would predict current positions.


## Our Questions:

1. Can we develop a simulation that accurately predicts planetary movements?
2. How do different numerical methods compare?
3. Are there limitations to the number of objects we can have in the simulation before it breaks down?


## Methods:
**System Class Methods:**
 - add_body
   - Adds bodies to our solar system in a discrete manner, allows free manipulation of what is and isn't included in the model
 - return_positions
   - Returns the positions of all the bodies added to the system.
 - return_bodies
   - This returns a pandas dataframe that displays all the data we have per body in the system.
   - Position, mass, velocity, diameter, gravity, and the index of the position
  - update_kinematics
    - pass
  - delete_body
    - Removes a body from the system.
  - interactive
    - This method creates a visualization of the data we have from return_bodies in two different forms:
      - 1d representation: This returns two graphs of a 1d representation of the magnitude of all the planets positions relative  to the sun, one graph returns the raw magnitude on a straight line, the other graph returns a logarithmic representation of the magnitude of the planets to create a more equally spaced representation to make parts of the original graph easier to read.
      - 3d representation: This returns a 3d representation of the positions of the planets relative to the sun given the data from return_bodies. X,Y,Z components are all represented here so you are able to visualize the relative 3d positions in relation to each other.
 - idealStepTimeValues
 - generate_SimulationData
 - generate_SimulationOrbitGraph
 - generate_SimulationAnimation


# Running Our Code:

**To run our code has a few distinct but simple steps.**

1. Make sure all the classes are properly implemented and defined.
2. Use the system.add_body to add planets to the solar system.
3. Use system.return_bodies to return a dataframe off all of the data added from step 2 to verify that data is accurate.
4. Use system.interactive(1)/system.interactive(3) to return interact-able 1d or 3d plots respectively to give a representation of your data from steps 2 and 3.
5. Animation stuff



# Conclusions
