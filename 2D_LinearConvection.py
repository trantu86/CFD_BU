##########################################################
# 2D Linear convection      du/dt + c*du/dx + c*du/dy = 0
#
##########################################################

from mpl_toolkits.mplot3d import axes3d
import numpy as np
import matplotlib.pyplot as plt
import time

nx = 20;    ny = 20;
nt = 50;    dt = 0.01;
c = 1;
dx = 2/(nx-1.0); dy = 2/(ny-1.0);

x = np.linspace(0,2,nx);
y = np.linspace(0,2,ny);
u = np.zeros((nx,ny), dtype=float);
uzero = np.zeros((nx,ny), dtype=float);

# initial conditions
for i in range(nx):
    for j in range(ny):
        if 0.5 <= x[i] <= 1 and 0.5 <= y[j] <= 1:
            u[i,j] = 2;
        else:
            u[i,j] = 1;

uzero[:] = u[:]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(x,y)
ax.plot_wireframe(X,Y,uzero)

for it in range(nt):
    un = u;
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            # FD in time BD in space
            u[i,j] = un[i,j] - c*dt/dx*(un[i,j] - un[i-1,j]) - c*dt/dy*(un[i,j]-un[i,j-1]);
        
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(x,y)
ax.plot_wireframe(X,Y,u)

plt.show()