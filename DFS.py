# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 21:44:45 2020

@author: David
"""
import numpy as np
import matplotlib.pyplot as plt
#init DFS 

g=9.81 #m/s
unit_x=np.array([1.,0.,0.])
unit_y=np.array([0.,1.,0.])
unit_z=np.array([0.,0.,1.])

class Disc:
    def __init__(self,mass,radius):
        self.mass = mass
        self.radius = radius
        self.pos_track=[]
        self.vel_track=[]

    def _get_force(self):
        F_g = -self.mass*g*unit_z
        return F_g
    
    def _get_acceleration(self):
        net_force = self._get_force()
        return net_force / self.mass

    def _update_pos_vel(self,dt):
        acceleration = self._get_acceleration()
        
        self.position=np.copy(self.position+dt*self.velocity)
        self.velocity=np.copy(self.velocity + dt*acceleration)

    def throw_disc(self,throw_v0,throw_p0,dt,t_thresh):
        t=0
        self.velocity = throw_v0
        self.position = throw_p0
        while t < t_thresh:
            #print(self.position)
            #print()
            self.pos_track.append(self.position)
            self.vel_track.append(self.velocity)
            self._update_pos_vel(dt)
            t+=dt
        #return np.array(pos_track),np.array(vel_track)
myDisc = Disc(100,10)

v0=np.array([0,0,100.])
p0=np.array([0,0,10.])
myDisc.throw_disc(v0,p0,.125,10)
zvals = np.array(myDisc.pos_track)[:,2]
plt.plot(zvals)
#print(myDisc.pos_track)
#plt.plot(pos_track[:,2])
#print(vel_track)