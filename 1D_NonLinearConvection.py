##########################################################
# 1D Non-Linear Convection      du/dt + u*du/dx = 0
# inviscid Burgers' equation
##########################################################

import numpy as np
import matplotlib.pyplot as plt
import time

### FD in time BD in space
def FDtimeBDspace(un,dt,dx,nx):
    u[:] = un[:]
    u[1:] = un[1:] - un[1:]*dt/dx*(un[1:] - un[0:-1])
#    for i in range(1,nx):        
#        u[i] = un[i] - un[i]*dt/dx*(un[i] - un[i-1])
        
    return u

### Lax-Friedrich scheme (1st order method) - Forward time centered scheme (FTCS) 
#   exact solution for sigma = 2.0 and shock wave of amplitude 1
def LFscheme(un,dt,dx,nx):
    u[:] = un[:]
    u[1:-1] = 1./2*(un[0:-2]+un[2:]) - 1./2*(un[0:-2]+un[2:])*dt/(2*dx)*(un[2:] - un[0:-2])
#    for i in range(1,nx-1):
#        u[i] = 1./2*(un[i-1]+un[i+1]) - 1./2*(un[i-1]+un[i+1])*dt/(2*dx)*(un[i+1] - un[i-1])

    return u

### Lax-Wendroff scheme (2nd order method)
def LWscheme(un,dt,dx,nx):
    u[:] = un[:]
    for i in range(1,nx-1):
        u[i] = un[i] - dt/4/dx*(un[i+1]**2-un[i-1]**2) + dt**2/4/dx**2*((un[i+1]+un[i])*(un[i+1]**2-un[i]**2)/2 - (un[i]+un[i-1])*(un[i]**2-un[i-1]**2)/2)

    return u
    
### MacCormack Scheme (2nd order method)
def MacCormack(un,dt,dx,nx):
    us = np.zeros(nx)
    u[:] = un[:]
    us[:] = un[:]
    for i in range(1,nx-1):
        us[i] = un[i] - dt/dx*(un[i+1]**2-un[i]**2)/2
        u[i] = 1./2*(un[i] + us[i] - dt/dx*(us[i]**2-us[i-1]**2)/2)
        
    return u

#------------------------------------------------------------------------------#
### parameters
xdim = 4;
nx = 41;                nt = 36;
sigma = 1.0;
dx = xdim/(nx-1.0);     dt = sigma*dx;

### arrays
x = np.linspace(0,xdim,nx)
u = np.zeros(nx)
un = np.zeros(nx)
uzero = np.zeros(nx)

# initial conditions
for i in range(nx):
    if 0.0 <= x[i] <= 2.0:
        u[i] = 1
    else:
        u[i] = 0

uzero[:] = u[:]

### plot initial condition
plt.plot(x,uzero,'b-')
line, = plt.plot(x,uzero,'k.-')         # Line2D object
plt.axis([-0.2, xdim+0.2, -0.2, 2.2])
plt.xlabel('x')
plt.ylabel('y')


### time evolution
for it in range(nt):
    un[:] = u[:]

#    u = FDtimeBDspace(un,dt,dx,nx)
#    u = LFscheme(un,dt,dx,nx)
#    u = LWscheme(un,dt,dx,nx)
    u = MacCormack(un,dt,dx,nx)
    
    line.set_ydata(u)       # Set the data np.array for y
    plt.draw()              # Redraw the current figure.
    plt.pause(0.1)


plt.hold(True)
