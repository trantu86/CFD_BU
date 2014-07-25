##########################################################
# Beam-Warming implicit method / Inciscid Burgers equation
#
##########################################################

from pylab import *
from numpy import *
from math import *
import matplotlib.pyplot as plt
#rcParans['font.family'] = 'serif'

### Parameters
domain = 4.
nx = 41
dx = domain/(nx-1)
T = 1.8
courant = 1.        # courant number
dt = courant*dx     # time step size
nt = int(T/dt)
Ee = 0.125          # damping coefficient
print('dx='+str(dx)+', dt='+str(dt))

### Initialisation: two sets of arrays, for undamped and damped cases
# for undamped BW:
x = arange(0,domain+dx,dx)
u1=zeros(nx); un1=zeros(nx);
a1=zeros(nx); b1=zeros(nx); c1=zeros(nx); d1=zeros(nx);
c12=zeros(nx); d12=zeros(nx); # to preserve diagonals for next time step

# for damped BW:
u2=zeros(nx); un2=zeros(nx);
a2=zeros(nx); b2=zeros(nx); c2=zeros(nx); d2=zeros(nx);
c22=zeros(nx); d22=zeros(nx);

### Initial conditions
for i in range(nx):
    if 0<=x[i] and x[i]<=2:
        u1[i]=1
        u2[i]=1

### time evolution
for it in range(nt):
    
# Solution without damping
    un1[:] = u1[:]
    a1[0] = -dt/dx/4*un1[0] # commenting this line will cause "drooping left"
    for i in range(1,nx):
        a1[i] = -dt/dx/4*un1[i-1]
    for i in range(nx):
        b1[i] = 1
    for i in range(nx-1):
        c1[i] = dt/dx/4*un1[i+1]
    for i in range(nx):
        d1[i] = un1[i]
        
    d1[0] = d1[0] - a1[0]*un1[0]  # take the left BC into account
    c12[0] = c1[0]/b1[0]
    d12[0] = d1[0]/b1[0]
    
    for i in range(1,nx-1):
        c12[i] = c1[i]/(b1[i]-c12[i-1]*a1[i])
        
    for i in range(1,nx):
        d12[i] = (d1[i]-d12[i-1]*a1[i])/(b1[i]-c12[i-1]*a1[i])
    
    u1[nx-1] = d12[nx-1]
    
    for i in range(nx-2,-1,-1):
        u1[i] = d12[i]-c12[i]*un1[i+1]
        
    # BC
    u1[0] = 1

# Solution with damping:
    un2[:] = u2[:]
    a2[0] = -dt/dx/4*un2[0]  # commenting this line will cause "drooping left"
    for i in range(1,nx):
        a2[i] = -dt/dx/4*un2[i-1]
    for i in range(nx):
        b2[i] = 1
    for i in range(nx-1):
        c2[i] = dt/dx/4*un2[i+1]
        
    d2[0] = un2[0]
    d2[1] = un2[1]
    # here we introduce damping
    for i in range(2,nx-2):
        d2[i] = un2[i] -Ee*(un2[i+2]-4*un2[i+1]+6*un2[i]-4*un2[i-1]+un2[i-2])
    d2[nx-2] = un2[nx-2]
    d2[nx-1] = un2[nx-1]
    
    d2[0] = d2[0] - a2[0]*un2[0]  # take the left BC into account
    c22[0] = c2[0]/b2[0]
    d22[0] = d2[0]/b2[0]
    
    for i in range(1,nx-1):
        c22[i] = c2[i]/(b2[i]-c22[i-1]*a2[i])
        
    for i in range(1,nx):
        d22[i] = (d2[i]-d22[i-1]*a2[i])/(b2[i]-c22[i-1]*a2[i])
        
    u2[nx-1] = d22[nx-1]
    
    for i in range(nx-2,-1,-1):
        u2[i] = d22[i] - c22[i]*un2[i+1]
        
    # BC
    u2[0] = 1


### Ploting
p1, = plt.plot(x,u1,'g.-')
p2, = plt.plot(x,u2,'r.-')
plt.axis([0, 4, -0.5, 2])
plt.xlabel('x')
plt.ylabel('u')
plt.legend([p1,p2],['undamped','damped with e = '+str(Ee)])
plt.show()