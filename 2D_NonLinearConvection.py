##########################################################
# 2D Non-Linear convection      du/dt + u*du/dx + v*du/dy = 0
#                               dv/dt + u*dv/dx + v*dv/dy = 0
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
v = np.zeros((nx,ny), dtype=float);

for i in range(nx):
    for j in range(ny):
        if 0.5 <= x[i] <= 1 and 0.5 <= y[j] <= 1:
            u[i,j] = 2;
            v[i,j] = 2;
        else:
            u[i,j] = 1;
            v[i,j] = 1;

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(x,y)
ax.plot_wireframe(X,Y,u)

# Finite Difference calculation
for it in range(nt):
    un = u;
    vn = v;
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            # FD in time BD in space
            u[i,j] = un[i,j] - un[i,j]*dt/dx*(un[i,j] - un[i-1,j]) - vn[i,j]*dt/dy*(un[i,j]-un[i,j-1]);
            v[i,j] = vn[i,j] - un[i,j]*dt/dx*(vn[i,j] - vn[i-1,j]) - vn[i,j]*dt/dy*(vn[i,j]-vn[i,j-1]);
        
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(x,y)
ax.plot_wireframe(X,Y,u)

plt.show()