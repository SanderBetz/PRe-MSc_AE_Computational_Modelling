#================================================================= 
#
# AE2220-II: Computational Modelling 
# Code for work session 1 - Preparation
#
#=================================================================
# This code provides a base for computing the advection 
# equation for Z on a rectangular domain
#
# lines 24-31:   Input parameters 
# lines 87-107:  Implementation of finite-difference scheme.
#                This is based on an explicit space march in y.
#                ****** indicates where information is to be added.
# 
#=================================================================

import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#=========================================================
# Input parameters
#=========================================================
nozc       = 0.0                # Nozzle centre
nozw       = 1.0                # Nozzle width
imax       = 80                 # Number of mesh points in x
nmax       = 80                 # Number of mesh points in y
k          = 0.05               # Artificial dissipation parameter

maxl       = 50                 # maximum grid lines on plots
stride     = 1                  # Point skip rate for suface plot


#=========================================================
# Load the mesh coordinates
#=========================================================
x1     = -2.5                    # Forward boundary position
x2     =  2.5                    # Rear boundary position
y1     =  0.                     # lower boundary position
y2     = 10.0                    # Upper boundary position
dx     = (x2-x1)/(imax-1)        # Mesh spacing in x
dy     = (y2-y1)/(nmax-1)        # Mesh spacing in y
x      = np.zeros((imax,nmax))   # Mesh x coordinates
y      = np.zeros((imax,nmax))   # Mesh y coordinates

# Mesh coordinates
for n in range(0, nmax):
  x[:,n] = np.linspace(x1, x2, imax)

for i in range(0, imax):
  y[i,:] = np.linspace(y1, y2, nmax)


#=========================================================
# Load the flow velocities 
#=========================================================
u  =np.zeros((imax,nmax))
v  =np.zeros((imax,nmax))

for n in range(0, nmax):
  for i in range(0, imax):
    rx=x[i,n]
    ry=y[i,n]-15
    r=math.sqrt(rx*rx+ry*ry)
    th=math.atan2(rx,ry)
    u[i,n]     = 4*math.sin(th)/r
    v[i,n]     = 1.0 + 4*math.cos(th)/r
    

#=========================================================
# Initialise the mass fraction distribution to zero
# then set the inflow condition for Z
#=========================================================
Z = np.zeros((imax,nmax))

xnozLeft =nozc-0.5*nozw
xnozRight=nozc+0.5*nozw
for i in range(0, imax):
  if (x[i,0] >xnozLeft) & (x[i,0] <xnozRight) :
    delx=x[i,0]-nozc
    Z[i,0] = 1*math.exp(-10*delx*delx)


#=========================================================
# March explicitly in y, solving for 
# the local value of the mass fraction Z
#=========================================================

for n in range(0, nmax-1):   # March upwards from y=0 to y=10
 
   # Update left boundary (node i = 0)  ******
   # Since u<0 on the left boundary a (first-order) numerical condition
   # should be implemented. Finish the line and uncomment it
#  Z[0,n+1] = Z[0,n] - \
 
  
   # Update interior (nodes i=1 to imax-2) 
   # Here we use a second-order central approximation for du/dx 
   # with artificial viscosity.
   for i in range(1, imax-1):
     Z[i,n+1] = Z[i,n] - (Z[i+1,n]-Z[i-1,n])*dy*u[i,n]/(v[i,n]*2.*dx) + \
                       k*(Z[i+1,n]-2*Z[i,n]+Z[i-1,n]) 


   # Update right boundary (node i = imax-1)  ******
   # Since u>0 on the right boundary a (first-order) numerical condition
   # should be implemented. Finish the line and uncomment it.
#  Z[imax-1,n+1] =  Z[imax-1,n] - \



#=========================================================
#  Plot results
#=========================================================
fig = plt.figure(figsize=(8,8))

ax1 = plt.subplot2grid((2,2), (0,0), colspan=3, rowspan=2, projection='3d')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('Z')
ax1.plot_surface(x, y, Z, shade=True, rstride=stride,cstride=stride,
cmap=plt.cm.gnuplot, linewidth=0, antialiased=True)

ax1.view_init(30, -120)
#plt.savefig('flame.png',dpi=250)
plt.show()

print('done')
