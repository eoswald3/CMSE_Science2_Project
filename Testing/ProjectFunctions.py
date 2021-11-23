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
    Calculate the acceleration on an object due to n bodies aroun it
    
    objIndex: index of the object the force is acting on {int}
    vectObjPos: vectorized list of object positions {np.array of shape(n,m)}
        n: spacial dimensions of the model
        m: number of objects
    objMasses: list of object masses {np.array of shape(m)}
    
    return a vector containing the force acting on the object
    '''
    G = 6.67408E-11
    
    mass_WO_Obj = np.delete(objMasses, objIndex) #Mass of the objects not including the origin object
    pos_WO_Obj = np.delete(vectObjPos, objIndex, axis = 1) #Position of the objects not including the origin object
    pos_Obj = vectObjPos[:,objIndex] #Position of the origin object
    pos_Obj = np.reshape(pos_Obj, (len(pos_Obj),1)) #Reshape of a vector

    forceOnObject = np.sum(G * objMasses[objIndex] * mass_WO_Obj * (pos_WO_Obj - pos_Obj) * (1/(eucDist(pos_Obj, pos_WO_Obj))**3), axis=0)
    
    return(forceOnObject)
    
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