#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 13:16:01 2021

@author: steven
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px #import for 3d interactive plotting
import ProjectFunctions as pf

class system():
    
    def __init__(self, name, mass, diameter, pos=[0,0,0], gravity=0): #diameter in kilometers. 
        
        '''
        Initializes the class.
        
        Takes inputs for the central body of the system.
            Central body at pos (0,0,0)
        
        Creates Array and DataFrame of Info.
            
        '''
        
        #initializing the body dataframe
        self.bodies = pd.DataFrame(columns=['Initial Position (m)','Initial Velocity (m/s)', 'Mass (kg)', 'Diameter (m)', 
                                            'Gravity (m/s^2)', 'Position Index'])
        
        #numpy array of shape dimension 3, number of objects
        self.positions = np.array([pos]).T
        self.velocities = np.array([pos]).T
        self.masses = np.array(mass)
        self.names = np.array(name)
        
        self.central_pos = pos #setting the position of the central body
        self.central_name = name #setting the name of the central body
        self.central_mass = mass #setting the mass of the central body
        
        #adding the parameters of the central body to the dataframe
        self.bodies.loc[name] = {'Initial Position (m)': pos, 'Mass (kg)': mass, 'Initial Velocity (m/s)': [0,0,0], 
                                 'Diameter (m)': diameter, 'Position Index': 0} 
    
    
    def add_body(self, name, mass, pos, velocity, diameter, gravity=0):
        
        '''
        Adds an individual body into the system.
        '''
        
        #find the index of the new row
        index = len(self.bodies) 
        
        #create dict of data for new row
        new_row = {'Initial Position (m)': pos, 'Mass (kg)': mass, 'Initial Velocity (m/s)': velocity, 
                   'Diameter (m)': diameter, 'Gravity (m/s^2)': gravity, 'Position Index': index}
        
        #create array of positions of shape dimension 3 by 1
        pos = np.array([pos]).T 
        vel = np.array([velocity]).T
        
        #insert new row in dataframe
        self.bodies.loc[name] = new_row 
        
        #append individual position,velocity,and mass information to the respective full array
        self.positions = np.append(self.positions, pos, axis=1) 
        self.velocities = np.append(self.velocities, vel, axis=1)
        self.masses = np.append(self.masses, mass)
        self.names = np.append(self.names, name)
        
    def return_positions(self):
        
        '''
        Return fuction for positions
        '''
    
        return self.positions
    
    def return_bodies(self):
        
        '''
        Return fuction for the bodies DataFrame.
        '''
        
        return self.bodies
    
    def delete_body(self, body):
        
        '''
        Removes a single body from the system
        
        Body: The name of the body you wish to remove. Dtype = string.
        '''
        
        #pulls body's position index in the array from dataframe
        index = self.bodies.loc[body,'Position Index'] 
        
        #removes the body from the vectorized data arrays
        self.positions = np.delete(self.positions, index, 1) 
        self.velocities = np.delete(self.velocities, index, 1)
        self.masses = np.delete(self.masses, index)
        self.names = np.delete(self.names, index)
        
        #removes the body from the dataframe
        self.bodies = self.bodies.drop(body) 
        
        #resets the index position column in the DF to be in line with positions array.
        for row in range(len(self.bodies)):
            if self.bodies.iloc[row, 5] > index: 
                self.bodies.iloc[row, 5] -= 1
                
    
    def interactive(self, dimension):
        
        '''
        This method provides two different interactive looks at the positions at the most recent update.
        Dimension: The dimension in which they want to view the system. Options: 1, or 3. Input as an integer.
        
        This method can only be run after generate_SimulationData is run, because of calling positional values.
        '''
        
        #insert update_update kinematics function
        
        #calling and saving the diameter of each body to use later.
        diameter = self.bodies['Diameter (m)'].to_list()
        
        #sets a color map so each body has a different color based on their initial velocity
        c = self.velocities[0]
        
        
        '''
        This will return the three dimensional interactive graph.
        
        The if statments use the argument entered into the method.
        '''
        
        #If three dimensions are requested
        if dimension == 3:

            #initializing x, y, and z values
            x = self.sim_positions[len(self.sim_positions)-1, 0] 
            y = self.sim_positions[len(self.sim_positions)-1, 1]
            z = self.sim_positions[len(self.sim_positions)-1, 2]

            #create the 3d plot
            scatter3d = px.scatter_3d(x=x, y=y, z=z, color=c, text=self.bodies.index.to_list(), template='plotly_dark', 
                                      hover_name=self.bodies.index.to_list(), size = diameter, title = 'Current Positions and Data')
            
            #update the background, and aesthetics for a cleaner image.
            scatter3d.update_scenes(xaxis_showgrid=False,yaxis_showgrid=False, zaxis_showgrid=False, xaxis_showticklabels=False, 
                                    yaxis_showticklabels=False, zaxis_showticklabels=False, camera_projection_type='orthographic',
                                    xaxis_zeroline=False, yaxis_zeroline=False, zaxis_zeroline=False)
            
            #removes the colorbar from the image
            scatter3d.update_coloraxes(showscale=False)
            
            #makes the hoverdata display next to the body
            scatter3d.update_layout(hovermode='x')
            
            #returns the 3d scatter plot
            return scatter3d.show()
        
        
        '''
        This will return the one dimensional interactive graphs.
        
        It returns two graphs of the magnitude of the distance of each body from the central body.
         - Graph 1: The real scale of distance.
         - Graph 2: The logarithmic scale of distance.
        
        The if statments use the argument entered into the method.
        '''
        
        # if one dimension is requested
        if dimension == 1:
            
            #initializing the list of magnitudes
            mag = []
            
            #Looping over the planets in the system by index of the planet
            for i in range(0,len(self.bodies)):
                planet_mag = np.linalg.norm(self.positions[:,i]) # np.linalg.norm takes in an multidimensional array and returns the magnitude of that vector.
                mag.append(planet_mag) # adds the magnitudes to a new list.
        
            # graphing the magnitudes on a 1d plane being labeled the same way as the 3d interactive graph. Y axis is simply zeros.
            fig = px.scatter(x=mag, y=np.zeros(11),hover_name=self.bodies.index.to_list(), color=c, 
                            size=diameter, template='plotly_dark', title = 'Magnitude of Distance Normal') 

            # same as before but we want our x axis to be logarithmic so some bodies are easier to view.
            logfig = px.scatter(x=mag, y=np.zeros(11),hover_name=self.bodies.index.to_list(),log_x=True, 
                                color=c, size=diameter, template='plotly_dark', title = 'Magnitude of Distance Logarithmic')
            
            # returns both figures
            return fig.show(), logfig.show()
        
        # if an integer other than 1 or 3 is inputed, it will return this error message.
        else:
            
            print('Uh Oh! That dimension is not correct! Please input 1 or 3.') # will print the two argument options.
            
    def calculate_residuals(self, correct_pos):
        
            '''
            This method would only be used when using actual planets/bodies.
            
            correct_pos: the correct positions of the planets based on actual data.
            
            returns an array of residuals that are the difference
            between the correct plantery positions and the predicted ones.
            '''
            
            for i in range(len(self.sim_positions)):
                val = np.linalg.norm(self.sim_positions[i])
            residuals = correct_pos - val
            return residuals
            
    
    def idealStepTimeValues(self, numYears):
        iterTime = 365*numYears
        h = 84000/numYears
        t = np.arange(0,iterTime,1)
        
        print('To model', numYears, 'years with a step-size (h) value of 1, set h and time parameters to the following.')
        print('h:',h)
        print('time:',iterTime)
    
    def generate_SimulationData(self, h, time=1, method=0):
        '''
        h: number of iterations in the simulation
        time: number of _____ to run the simulation for !!!
        method: method for generating data
                0: Euler's Method
                1: Huen's Method
                2: Verlat's Method
                etc...
                
        return a vectorized array of positions, velocities, and accelerations for each body over 
               the course of the simulation using the chosen method
        '''
        time = int(time)
        
        m_aU = 1/(1.494E11) #convert from meters to aU
        mS_aUYr = 1/4744 #convert from meters/sec to aU/yr
        earthMass = 1/(5.972E24) #Convert from kg to Earth masses
        
        if method == 0:
            self.sim_positions, self.sim_velocities, self.sim_accelerations = pf.kinematic_Euler(h,time,self.positions,self.velocities,self.masses)
 
        elif method == 1:
            self.sim_positions, self.sim_velocities, self.sim_accelerations = pf.kinematic_Huens(h,time,self.positions,self.velocities,self.masses)

        elif method == 2:
            h2 = h/3.145e7 #Convert from seconds to years
            
            self.sim_positions, self.sim_velocities, self.sim_accelerations = pf.kinematic_Verlat(h2,time,self.positions*m_aU,self.velocities*mS_aUYr,self.masses*earthMass)
            self.sim_positions = self.sim_positions/m_aU #Convert back to mks units
            self.sim_velocities = self.sim_velocities/mS_aUYr #Convert back to mks units
            self.sim_accelerations = self.sim_accelerations/earthMass #Convert back to mks units
        
        else:
            raise ValueError(method,' Is not a valid method identifier, please input a valid identifier')
            
    def generate_SimulationOrbitGraph(self, path, prefix, title, axesTF = True, second_View = False):
        '''
        generate a graph containing the plots of the objects path over the course of the simulation
        
        path: relative path to were the animation frames will be stored (end with /)
        prefix: generated image filename prefix
        axesTF: conditional for whether to have the axes elements on or off
        second_View: conditional for whether to plot a second view, orthogonal to the first
        '''
        plt.style.use('default') #Set style of the plot
        
        #Check which type of image is to be created
        if second_View == True:
            fig = plt.figure(figsize = (10,15)) #Set figure size
            ax = fig.add_subplot(211, projection='3d')
            
        else:
            fig = plt.figure(figsize = (10,10)) #Set figure size
            ax = fig.add_subplot(111, projection='3d')
            
        names = self.bodies.index.to_list() #Get names of objects to be ploted 
        
        for i in range(np.shape(self.sim_positions)[2]):
            xPlot = self.sim_positions[:,0,i] # Grab x position data of object i
            yPlot = self.sim_positions[:,1,i] # Grab y position data of object i
            zPlot = self.sim_positions[:,2,i] # Grab z position data of object i

            ax.plot(xPlot, yPlot, zPlot, label = names[i]) #Plot the full orbit of each object
            ax.scatter(xPlot[0], yPlot[0], zPlot[0], c = 'green') #Indicate initial position
            ax.scatter(xPlot[-1], yPlot[-1], zPlot[-1], c = 'red') #Indicate final position

        ax.view_init(elev=90, azim=90) #Change view parameters (spin along the azimuth)
        
        ax.set_xlabel('x', fontsize = 15) #Assign label
        ax.set_ylabel('y', fontsize = 15) #Assign label
        ax.set_zlabel('z', fontsize = 15) #Assign label
        
        #Plot second view if conditional is True
        if second_View == True:
            ax1 = fig.add_subplot(212, projection='3d')

            for i in range(np.shape(self.sim_positions)[2]):
                xPlot = self.sim_positions[:,0,i] # Grab x position data of object i
                yPlot = self.sim_positions[:,1,i] # Grab y position data of object i
                zPlot = self.sim_positions[:,2,i] # Grab z position data of object i

                ax1.plot(xPlot, yPlot, zPlot, label = names[i]) #Plot the full orbit of each object
                ax1.scatter(xPlot[0], yPlot[0], zPlot[0], c = 'green') #Indicate initial position
                ax1.scatter(xPlot[-1], yPlot[-1], zPlot[-1], c = 'red') #Indicate final position
                
            ax1.view_init(elev=0, azim=90) #Change view parameters (spin along the azimuth)
            
            ax1.set_xlabel('x', fontsize = 15) #Assign label
            ax1.set_ylabel('y', fontsize = 15) #Assign label
            ax1.set_zlabel('z', fontsize = 15) #Assign label
            
            plt.subplots_adjust(hspace=-0.3)

        #Set axes state
        ax.axis('on')
        if axesTF == False:
            ax.axis('off')

        #systemScale = 10e13

        #ax.set_xlim3d(-systemScale,systemScale)
        #ax.set_ylim3d(-systemScale,systemScale)
        #ax.set_zlim3d(-systemScale,systemScale)

        ax.set_title(title, y=0.95, fontsize=15)
        
        ax.legend()

        #Saving the graph
        fig.savefig(path + prefix + '.png',transparent=True)
        fig.clear()
        plt.close(fig)
        
    def generate_SimulationAnimation(self, path, prefix, title, scale = 1e11, axesTF = True):
        '''
        run the simulation with given input parameters and save the output images
        
        path: relative path to were the animation frames will be stored (end with /)
        prefix: generated image filename prefix
        title:
        scale: scale of the plot
        axesTF: conditional for whether to have the axes elements on or off
        
        #ffmpeg -framerate 10 -i Test_3DScatter_s40_cross_%04d.jpeg  Test_3DScatter_s40_3d_fr_cross_Movie.mp4
        '''
        
        #Use generated data to produce an animation for the same data
        xPlot, yPlot, zPlot = self.sim_positions[:,0], self.sim_positions[:,1], self.sim_positions[:,2]
        
        num_iterations = len(xPlot)
        
        diameter = self.bodies['Diameter (m)'].to_list()
        names = self.bodies.index.to_list()
    
        plt.style.use('dark_background')
        
        #Generate, output, and clear the figure for each frame
        count = 0 #Used to update the image filenames
        for i in range(0,num_iterations,1):
            fig = plt.figure(figsize = (10,10))
            ax = fig.add_subplot(111, projection='3d')
            
            diameterShow = (np.round(diameter/np.min(diameter))*5)
            diameterShow[0] = np.max(diameterShow)*1.5
            
            #tempAx = ax.scatter(xPlot[i], yPlot[i], zPlot[i], c = c, s = np.round(diameter/np.min(diameter))*4, edgecolors = 'black', alpha = 1, label = c)
            
            for k in range(len(xPlot[i])):
                ax.scatter(xPlot[i,k], yPlot[i,k], zPlot[i,k], s = diameterShow[k], edgecolors = 'white', alpha = 1, label = names[k])
            
            ax.view_init(elev=45, azim=-90) #Change view parameters (spin along the azimuth)
            
            #Set axes state
            ax.axis('on')
            if axesTF == False:
                ax.axis('off')
            
            ax.set_xlabel('x', fontsize = 15) #Assign label
            ax.set_ylabel('y', fontsize = 15) #Assign label
            ax.set_zlabel('z', fontsize = 15) #Assign label
            
            ax.set_xlim3d(-scale,scale)
            ax.set_ylim3d(-scale,scale)
            ax.set_zlim3d(-scale,scale)
            
            ax.set_title(title, y=0.95, fontsize=15)
            
            ax.legend(borderpad=1, labelspacing=1)
            
            #Saving of each frame
            i_str = str(count)
            suffix = i_str.rjust(4,'0')
            fig.savefig(path + prefix + suffix + '.jpeg')
            fig.clear()
            plt.close(fig)
            count += 1