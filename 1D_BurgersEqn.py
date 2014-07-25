##########################################################
# 1D Burgers' equation      du/dt + u*du/dx = vis*d^2u/dx^2
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
ua = np.zeros(nx);
ipl = [];
iml = [];

for i in range(nx):
    ipl.append(i+1);
    iml.append(i-1);
    x[i] = (i)*dx;
    
ipl[nx-1] = 0;
iml[0] = nx-1;

for i in range(nx):
    phi = np.exp(-x[i]**2/(4*vis)) + np.exp(-(x[i]-2*np.pi)**2/(4*vis));
    dphi = -0.5/vis*x[i]*np.exp(-x[i]**2/(4*vis)) - 0.5/vis*(x[i]-2*np.pi)*np.exp(-(x[i]-2*np.pi)**2/(4*vis));
    u[i] = -2*vis*dphi/phi + 4;

plt.plot(x,u);
#plt.show();

for it in range(nt):
    t = it*dt;
    for i in range(nx):
        phi = np.exp(-(x[i]-4*t)**2/(4*vis*(t+1))) * np.exp(-(x[i]-4*t-2*np.pi)**2/(4*vis*(t+1)));
        dphi = -0.5/(vis*(t+1))*(x[i]-4*t)*np.exp(-(x[i]-4*t)**2/(4*vis*(t+1))) - 0.5/(vis*(t+1))*(x[i]-4*t-2*np.pi)*np.exp(-(x[i]-4*t-2*np.pi)**2/(4*vis/(t+1)));
        ua[i] = -2*vis*dphi/phi + 4;
        
    un = u;
    for i in range(nx):
        u[i] = un[i] - un[i]*dt/dx*(un[i] - un[iml[i]]) + vis*dt/dx**2*(un[ipl[i]]-2*un[i]+un[iml[i]]);
    
#plt.plot(x,u);
#plt.plot(x,ua);
plt.show();