# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 21:44:45 2020

@author: David
"""
import numpy as np
import matplotlib.pyplot as plt
#init DFS 

class Disc:
    def __init__(self,mass,radius):
        self.mass = mass
        self.radius = radius
        self.velocity=np.array([0,0,0]) 
        self.position=np.array([0,0,0])
    
        

testDisc=Disc(10,5)
print(testDisc.velocity)