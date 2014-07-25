##########################################################
# 1D Diffusion      du/dt = vis*d^2u/dx^2
#
##########################################################

import numpy as np
import matplotlib.pyplot as plt
import time

nx = 20;    nt = 50;
dt = 0.01;  vis = 0.1;
dx = 2/(nx-1.0);

x = np.linspace(0,2,nx);
u = np.zeros(nx);
un = np.zeros(nx);
uzero = np.zeros(nx)

# initial conditions
for i in range(nx):
    if 0.5 <= x[i] <= 1:
        u[i] = 2;
    else:
        u[i] = 1;
        
uzero[:] = u[:]

line, = plt.plot(x,uzero,'k.-')         # Line2D object
#plt.axis([0, 2, -2, 2])
plt.xlabel('x')
plt.ylabel('y')

for it in range(nt):
    un[:] = u[:];
    for i in range(2,nx-1):
        # FD in time CD in space
        u[i] = un[i] + vis*dt/dx/dx*(un[i+1]-2*un[i]+un[i-1]);
        
    line.set_ydata(u)       # Set the data np.array for y
    plt.draw()              # Redraw the current figure.
    plt.pause(0.1)


plt.hold(True)
plt.plot(x,uzero,'b-')