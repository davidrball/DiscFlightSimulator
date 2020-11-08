# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 21:44:45 2020

@author: David
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#init DFS 

g=9.81 #m/s
rho_air = 1.22 #kg/m^3
C_d = 1. #drag coefficient (standin val for now)
C_l = 1. #lift coeff
unit_x=np.array([1.,0.,0.])
unit_y=np.array([0.,1.,0.])
unit_z=np.array([0.,0.,1.])
ground_z=0

#generate rotation matrics for given angles
#not the correct way to handle this, just use spherical coordinates theta and phi
def Rx(theta):
    return np.array([[1,0,0],[0,np.cos(theta),-np.sin(theta)],[0,np.sin(theta),np.cos(theta)]])
def Ry(theta):
    return np.array([[np.cos(theta),0,np.sin(theta)],[0,1,0],[-np.sin(theta),0,np.cos(theta)]])
yrot=45*np.pi/180.
xrot=45*np.pi/180.
Mxrot = Rx(xrot)
Myrot = Ry(yrot)
#rotated_vec =  np.dot(Myrot,np.dot(Mxrot,unit_z))
#totrot = np.dot(Myrot,Mxrot)



rotated_vec = np.dot(Mxrot,unit_z) #apply pitch first
rotated_vec = np.dot(np.dot(Mxrot,Myrot),rotated_vec)

print(rotated_vec)




class Disc:
    def __init__(self,mass,radius,height):
        self.mass = mass
        self.radius = radius
        self.height=height
        self.area=2*radius*height #assuming edge-on
        self.I=mass*radius**2 #assumes all mass is on rim (probably close)
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

    def throw_disc(self,throw_v0,throw_p0,throw_pitch0,throw_roll0,throw_omega0,dt,t_thresh):
        t=0
        self.velocity = throw_v0
        self.position = throw_p0
        self.L=throw_omega0*self.I
        self.pitch = throw_pitch0
        self.roll = throw_roll0
        
        
        #Mxrot=Rx(throw_rot0[0])
        #Myrot=Ry(throw_rot0[1])
        #spin_unitvec = np.dot(Myrot,np.dot(Mxrot,unit_z))
        
        #self.Lvec=L*()
        #self.Lvec=self.L*
        #Lvec_pitch = np.sin(throw_rot0[0])*unit_y
        #Lvec_roll = np.sin(throw_rot0[1])*unit_x
        #Lvec_z = (1-np.sqrt(Lvec_pitch**2+Lvec_roll**2))*unit_z
        
        
        while self.position[2] >ground_z:
            #print(self.position)
            #print()
            self.pos_track.append(list(self.position))
            self.vel_track.append(list(self.velocity))
            self._update_pos_vel(dt)
            t+=dt
        #return np.array(pos_track),np.array(vel_track)
        self.pos_track=np.array(self.pos_track) #just for easier manipulation at end
        self.vel_track=np.array(self.pos_track)
        
        
        
    def plot_trajectory(self):
        fig=plt.figure()
        ax=fig.add_subplot(111,projection='3d')
        ax.plot(np.array(self.pos_track)[:,0],np.array(self.pos_track)[:,1],np.array(self.pos_track)[:,2])
        ax.scatter(self.pos_track[0,0],self.pos_track[0,1],self.pos_track[0,2],color="C1")
# =============================================================================
# 
# mass=.17 #kg
# radius=.05 #cm
# height=0.01
# 
# tfinal=50
# dt=.125
# myDisc = Disc(mass,radius,height)
# 
# v0=np.array([0,70.,0.])
# p0=np.array([0,0,1.])
# rot0=np.array([0.,0.])*np.pi/180. #rotation about x (pitch) and about y (roll)
# omega0=50 #rad/s, sign indicates whether it's in +/- z direction
# myDisc.throw_disc(v0,p0,rot0,omega0,dt,tfinal)
# myDisc.plot_trajectory()
# 
# =============================================================================
