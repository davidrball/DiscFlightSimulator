# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 21:44:45 2020

@author: David
"""
import numpy as np
import matplotlib.pyplot as plt
#init DFS 

g=9.81 #m/s
rho_air = 1.22 #kg/m^3
C_d = 1. #drag coefficient (standin val for now)
C_l = 1. #lift coeff
unit_x=np.array([1.,0.,0.])
unit_y=np.array([0.,1.,0.])
unit_z=np.array([0.,0.,1.])
ground_z=0

class Disc:
    def __init__(self,mass,radius,height):
        self.mass = mass
        self.radius = radius
        self.height=height
        self.area=2*radius*height #assuming edge-on
        self.pos_track=[]
        self.vel_track=[]

    


    def _get_force(self):
        F_g = -self.mass*g*unit_z #gravity
        F_drag = -0.5*np.dot(self.velocity,self.velocity)*C_d*rho_air*self.area*self.velocity/np.sqrt(np.dot(self.velocity,self.velocity))
        F_lift = .5*np.dot(self.velocity,self.velocity)*C_l*rho_air*self.area*unit_z #for now just straight up
        
        F_net=F_g+F_drag+F_lift
        return F_net
    
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
        while self.position[2] >ground_z:
            #print(self.position)
            #print()
            self.pos_track.append(self.position)
            self.vel_track.append(self.velocity)
            self._update_pos_vel(dt)
            t+=dt
        #return np.array(pos_track),np.array(vel_track)

mass=.17 #kg
radius=.05 #cm
height=0.01

tfinal=50
dt=.125
myDisc = Disc(mass,radius,height)

v0=np.array([0,30.,5])
p0=np.array([0,0,10.])
myDisc.throw_disc(v0,p0,dt,tfinal)
zvals = np.array(myDisc.pos_track)[:,2]
yvals = np.array(myDisc.pos_track)[:,1]
plt.plot(zvals)
plt.plot(yvals)
#print(myDisc.pos_track)
#plt.plot(pos_track[:,2])
#print(vel_track)