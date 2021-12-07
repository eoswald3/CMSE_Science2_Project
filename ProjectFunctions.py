# -*- coding: utf-8 -*-
"""
Author: Steven VanCamp
Project: CMSE 201 Project
Co-Authors: Patrick, Ethan, Emma 
"""
import numpy as np

def eucDist(origin, vectorizedGrid):
    '''
    Calculate the euclidian distance between a point (origin) and another point(s)
    
    origin: origin points np.array of shape{(n,1)}
        n: spacial dimensions of the model
    vectorizedGrid: vectorized grid {np.array of shape(n,m^n)} ????
        m: grid resolution ????
        
    return an array containing the distance to each point on the grid from the origin
    '''
    dist = np.sqrt(np.sum((origin-vectorizedGrid)**2, axis=0))

    return(dist)

def fofGrav(objIndex, vectObjPos, objMasses):
    '''
    Calculate the acceleration on an object due to n bodies around it
    
    objIndex: index of the object the force is acting on {int}
    vectObjPos: vectorized list of object positions {np.array of shape(n,m)}
        n: spacial dimensions of the model
        m: number of objects
    objMasses: list of object masses {numpy.array of shape (m)}
    
    return a vector containing the force acting on the object
    '''
    G = 6.67408E-11
    
    mass_WO_Obj = np.delete(objMasses, objIndex) #Mass of the objects not including the origin object
    pos_WO_Obj = np.delete(vectObjPos, objIndex, axis = 1) #Position of the objects not including the origin object
    pos_Obj = vectObjPos[:,objIndex] #Position of the origin object
    pos_Obj = np.reshape(pos_Obj, (len(pos_Obj),1)) #Reshape of a vector

    forceOnObject = np.sum(G * objMasses[objIndex] * mass_WO_Obj * (pos_WO_Obj - pos_Obj) * (1/(eucDist(pos_Obj, pos_WO_Obj))**3), axis=1)
    
    return(forceOnObject, forceOnObject/objMasses[objIndex])

def acceleration(objIndex, vectObjPos, objMasses, G = 4*(np.pi**2/333000)):
    '''
    Alternate method for calculating acceleration on an object due to n bodies around it
    
    objIndex: index of the object the force is acting on {int}
    vectObjPos: vectorized list of object positions {np.array of shape(n,m)}
        n: spacial dimensions of the model
        m: number of objects
    objMasses: list of object masses {np.array of shape(m)}
    
    return a vector containing the accelerations on the object of {numpy.array of shape (3,1)}
    '''
    pos_WO_Obj = np.delete(vectObjPos, objIndex, axis = 1) #Position of the objects not including the origin object
    pos_Obj = vectObjPos[:,objIndex] #Position of the origin object
    pos_Obj = np.reshape(pos_Obj, (len(pos_Obj),1)) #Reshape of a vector
    
    mass_WO_Obj = np.delete(objMasses, objIndex) #Masses of all objects, not including the origin object
    
    radialDist = np.linalg.norm(pos_Obj-pos_WO_Obj,axis=0) #Calculate the radial distance between the origin object and all orther objects
    
    a = -G*(np.sum(mass_WO_Obj*(pos_Obj-pos_WO_Obj)*(1/np.abs(radialDist)**3),axis=1)) #Calculate the acceleration on the origin object

    a = np.reshape(a,np.shape(pos_Obj)) #Reshape the acceleration array into the correct form
    
    return(a)
    
def vectorizeGrid(grid): #Used for np.mgrid or a combined np.meshgrid
    '''
    Vectorize and decompose the components of the grid
    
    grid: np.array of shape(n,m...m_n)
        n: spacial dimensions of the model
        m: grid resolution
        
    return an array containing the vectorized components of grid shape(n,m^n) ????
    '''
    dimension = grid.shape[0] #Get n
    resolution = np.prod(grid.shape[1:]) #Get resolution by taking the product of the m values 
    
    vectorizedGrid = np.reshape(grid[:],(dimension,resolution)) #Get first vectorized component

    return(vectorizedGrid)

def kinematic_Euler(h,iterations,pos_Array,vel_Array,mass_Array):
    '''
    h: stepsize
    iterations: number of iterations
    pos_Array: array of initial positions (x,y,z) [m] of shape ()
    vel_Array: array of initial velocites [m/s] of shape ()
    mass_Array: array of masses [kg] of shape ()
    
    return a vectorized list of positions, velocities, and accelerations for each body over the course of the simulation using Euler's method
    '''
    #Set up array of times to iterate through
    t = np.arange(0,iterations,1)
    
    #Get initial acceleration
    acc_Array = fofGrav(0,pos_Array,mass_Array)[1]

    for i in range(pos_Array.shape[1]-1):
        acc_Array = np.vstack((acc_Array, fofGrav(i+1,pos_Array,mass_Array)[1]))
    
    acc_Array = acc_Array.T
    
    #Set up arrays for value storing/output
    nObjects = pos_Array.shape

    v = np.zeros((iterations,3,nObjects[1]))
    r = np.zeros((iterations,3,nObjects[1]))
    a = np.zeros((iterations,3,nObjects[1]))
    
    #Set storage arrays 0 position to initial values
    v[0] = vel_Array
    r[0] = pos_Array
    a[0] = acc_Array
    
    #Calculate position, velocity, and acceleration values
    for i in range(0,len(t)-1):
        v[i+1] = v[i] + h * a[i]
        r[i+1] = r[i] + h * v[i]
        
        #Get accelerations
        acc_Array = fofGrav(0,r[i],mass_Array)[1]
    
        #Calculate the acceleration of the remaining objects
        for k in range(pos_Array.shape[1]-1):
            acc_Array = np.vstack((acc_Array, fofGrav(k+1,r[i],mass_Array)[1]))
    
        a[i+1] = acc_Array.T
        
    return(r,v,a)

def kinematic_Huens(h,iterations,pos_Array,vel_Array,mass_Array):
    '''
    h: stepsize
    iterations: number of iterations simulation
    pos_Array: array of initial positions (x,y,z) [m] of shape ()
    vel_Array: array of initial velocites [m/s] of shape ()
    mass_Array: array of masses [kg] of shape ()
    
    return a vectorized array of positions, velocities, and accelerations for each body over the course of the simulation using Huen's method
    '''
    #Set up array of times to iterate through
    t = np.arange(0,iterations,1)
    
    #Get initial acceleration
    acc_Array = fofGrav(0,pos_Array,mass_Array)[1]

    for i in range(pos_Array.shape[1]-1):
        acc_Array = np.vstack((acc_Array, fofGrav(i+1,pos_Array,mass_Array)[1]))
    
    acc_Array = acc_Array.T
    
    #Set up arrays for value storing/output
    nObjects = pos_Array.shape

    v = np.zeros((iterations,3,nObjects[1]))
    r = np.zeros((iterations,3,nObjects[1]))
    a = np.zeros((iterations,3,nObjects[1]))
    
    #Set storage arrays 0 position to initial values
    v[0] = vel_Array
    r[0] = pos_Array
    a[0] = acc_Array
    
    #Calculate position, velocity, and acceleration values
    for i in range(0,len(t)-1):
        v[i+1] = v[i] + h * a[i]
        r[i+1] = r[i] + (h/2)*(v[i]+v[i+1])
        
        #Get accelerations
        acc_Array = fofGrav(0,r[i],mass_Array)[1]
    
        #Calculate the acceleration of the remaining objects
        for k in range(pos_Array.shape[1]-1):
            acc_Array = np.vstack((acc_Array, fofGrav(k+1,r[i],mass_Array)[1]))
    
        a[i+1] = acc_Array.T
        
    return(r,v,a)

def kinematic_Verlat(h,iterations,pos_Array,vel_Array,mass_Array):
    '''
    h: stepsize
    iterations: number of iterations simulation
    pos_Array: array of initial positions (x,y,z) [m] of shape ()
    vel_Array: array of initial velocites [m/s] of shape ()
    mass_Array: array of masses [kg] of shape ()
    
    return a vectorized array of positions, velocities, and accelerations for each body over the course of the simulation using Huen's method
    '''
    #Set up array of times to iterate through
    t = np.arange(0,iterations,1)
    
    #Get initial acceleration
    acc_Array = acceleration(0,pos_Array,mass_Array)

    for i in range(pos_Array.shape[1]-1):
        acc_Array = np.hstack((acc_Array,acceleration(i+1,pos_Array,mass_Array)))
    
    #Set up arrays for value storing/output
    nObjects = pos_Array.shape

    v = np.zeros((iterations,3,nObjects[1]))
    r = np.zeros((iterations,3,nObjects[1]))
    a = np.zeros((iterations,3,nObjects[1]))
    
    #Set storage arrays 0 position to initial values
    v[0] = vel_Array
    r[0] = pos_Array
    a[0] = acc_Array
    
    #Calculate position, velocity, and acceleration values
    for i in range(0,len(t)-1):
        r[i+1] = r[i] + h*v[i] + ((h**2)/2) * (a[i])
        
        #Calculate the acceleration of the first object
        acc_Array = acceleration(0,r[i+1],mass_Array)
        
        #Calculate the acceleration of the remaining objects
        for k in range(pos_Array.shape[1]-1):
            acc_Array = np.hstack((acc_Array,acceleration(k+1,r[i+1],mass_Array)))
        
    a[i+1] = acc_Array
    
    v[i+1] = v[i] + h*a[i] + (h/2) * (a[i+1]-a[i])
    
    return(r,v,a)
    