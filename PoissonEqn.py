##########################################################
# Poisson equation      d^2p/dx^2 + d^2p/dy^2 = b
#
##########################################################

from mpl_toolkits.mplot3d import axes3d
import numpy as np
import matplotlib.pyplot as plt
import time

# parameters
nx = 20;    ny = 20;     nit = 1000;    
dx = 2./(nx-1); dy = 1./(ny-1);

# arrays
x = np.linspace(0,2,nx);
y = np.linspace(0,2,ny);
p = np.zeros((nx,ny), dtype=float);
pn = np.zeros((nx,ny), dtype=float);
b = np.zeros((nx,ny), dtype=float);

# initial conditions
b[nx/4,ny/4] = 100;
b[nx*3/4,ny*3/4] = -100;

# figure initialization
plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

X, Y = np.meshgrid(x,y)
wframe = ax.plot_wireframe(X,Y,p)

# time iteration
for iit in range(nit):
    pd = p;

    # space iteration
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            # 2nd order CD in space
            p[i,j] = ((pd[i+1,j] + pd[i-1,j])*dy**2 + (pd[i,j+1] + pd[i,j-1])*dx**2 - b[i,j]*dx**2*dy**2)/(dx**2 + dy**2)/2;
            
    p[1:nx-1,0] = p[1:nx-1,1]
    p[1:nx-1,ny-1] = p[1:nx-1,ny-2]
    
    # Draw new frame
    if np.floor((iit+1)/10) < np.floor((iit+2)/10): 
        # Remove old line collection before drawing
        oldcol = wframe
        if oldcol is not None:
            ax.collections.remove(oldcol)

        # Draw new wireframe
        wframe = ax.plot_wireframe(X,Y,p)
        plt.draw()
        plt.pause(0.1)