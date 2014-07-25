##########################################################
# 1D Linear convection      du/dt + c*du/dx = 0
#
##########################################################

import numpy as np
import matplotlib.pyplot as plt
import time

### FD in time BD in space scheme
def FDtimeBDspace(un,c,dt,dx,nx):
    u[:] = un[:]
    u[1:] = un[1:] - c*dt/dx*(un[1:] - un[0:-1])
#    for i in range(1,nx):        
#        u[i] = un[i] - c*dt/dx*(un[i] - un[i-1])    
    
    return u

### FD in time FD in space scheme
def FDtimeFDspace(un,c,dt,dx,nx):
    u[:] = un[:]
    u[0:-1] = un[0:-1] - c*dt/dx*(un[1:] - un[0:-1])
#    for i in range(0,nx-1):
#        u[i] = un[i] - c*dt/dx*(un[i+1] - un[i])

    return u

### FD in time CD in space scheme
def FDtimeCDspace(un,c,dt,dx,nx):
    u[:] = un[:]
    u[1:-1] = un[1:-1] - c*dt/(2*dx)*(un[2:] - un[0:-2])
 #   for i in range(1,nx-1):
 #       u[i] = un[i] - c*dt/(2*dx)*(un[i+1] - un[i-1])

    return u
    
### Lax-Friedrich scheme - Forward time centered scheme (FTCS)
#   exact solution for sigma = 1.0
def LFscheme(un,c,dt,dx,nx):
    u[:] = un[:]
    u[1:-1] = 1./2*(un[0:-2]+un[2:]) - c*dt/(2*dx)*(un[2:] - un[0:-2])
#    for i in range(1,nx-1):
#        u[i] = 1./2*(un[i-1]+un[i+1]) - c*dt/(2*dx)*(un[i+1] - un[i-1])

    return u

### Lax-Wendroff scheme
def LWscheme(un,c,dt,dx,nx):
    u[:] = un[:]
    u[1:-1] = u[1:-1] - 1./2*c*dt/dx*(un[2:] - un[0:-2]) + 1./2*(c*dt/dx)**2*(un[2:]-2*un[1:-1]+un[0:-2])
#    for i in range(1,nx-1):
#        u[i] = u[i] - 1./2*c*dt/dx*(un[i+1] - un[i-1]) + 1./2*(c*dt/dx)**2*(un[i+1]-2*un[i]+un[i-1])

    return u

### Leapfrog scheme
def Leapfrogscheme(unn,un,c,dt,dx,nx):
    u[:] = un[:]
    for i in range(1,nx-1):
        u[i] = unn[i] - c*dt/dx*(un[i+1]-un[i-1])
        
    return u
        
#------------------------------------------------------------------------------#
### parameters
nx = 41
nt = 40
c = 1
sigma = 1.0
#dx = 2/(nx-1.0)
dx = 0.05
dt = sigma*dx/c

### arrays
x = np.linspace(0,2,nx)
#x = np.arange(0,2+dx,dx)
u = np.zeros(nx)
un = np.zeros(nx)
unn = np.zeros(nx)
uzero = np.zeros(nx)

### initial conditions 
## square wave
for i in range(nx):
    if 0.5 <= x[i] <= 1:
        u[i] = 2
    else:
        u[i] = 1

## triangle wave
#for i in range(nx):
#    if 0.9 <= x[i] <= 1:
#        u[i] = 10*(x[i]-0.9)
#    if 1.0 <= x[i] <= 1.1:
#        u[i] = 10*(1.1-x[i])

## square wave at beginning        
#u[0] = 1

uzero[:] = u[:]

### plot initial conditions
plt.plot(x,uzero,'b--')
line, = plt.plot(x,uzero,'k.-')         # Line2D object
plt.axis([-0.2, 2.2, -2.0, 4.0])
plt.xlabel('x')
plt.ylabel('y')

### time evolution #######################
for it in range(nt):
    unn[:] = un[:]
    un[:] = u[:]

#    u = FDtimeBDspace(un,c,dt,dx,nx)            # FD in time BD in space scheme
#    u = FDtimeFDspace(un,c,dt,dx,nx)            # FD in time FD in space scheme
#    u = FDtimeCDspace(un,c,dt,dx,nx)            # FD in time CD in space scheme
    u = LFscheme(un,c,dt,dx,nx)                 # Lax-Friedrich scheme
#    u = LWscheme(un,c,dt,dx,nx)                 # Lax-Wendroff scheme
    # Leapfrog scheme
#    if it == 0:
#        u = FDtimeBDspace(un,c,dt,dx,nx)        # starting scheme
#    else:
#        u = Leapfrogscheme(unn,un,c,dt,dx,nx)   # Leapfrog scheme
    

    line.set_ydata(u)       # Set the data np.array for y
    plt.draw()              # Redraw the current figure.
    plt.pause(0.1)

#------------------------------------------------------------------------------#
### translation of initial condition
nx = 1000
x = np.linspace(0,2,nx)
ua = np.zeros(nx)
ddx = dt*nt*c

## square wave
for i in range(nx):
    if 0.5+ddx <= x[i] <= 1+ddx:
        ua[i] = 2
    else:
        ua[i] = 1

## triangle wave
#for i in range(nx):
#    if 0.9+ddx <= x[i] <= 1.0+ddx:
#        ua[i] = 10*(x[i]-0.9-ddx);
#    if 1.0+ddx <= x[i] <= 1.1+ddx:
#        ua[i] = 10*(1.1+ddx-x[i]);    
        
plt.plot(x,ua,'b--')
plt.pause(0.1)
