"""
SOAC Exercise 2: Molenkamp test
Perform Molenkamp simulation using three different numerical schemes:
    (1) Euler forward in time and upwind in space
    (2) "Lax-Wendroff" scheme
    (3) The spectral method

T. Bovenschen
E.A. Melman
@ Utrecht University
"""

# Import libraries
from IPython import get_ipython
get_ipython().magic('reset -sf')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import animation

# --------------------------------------------
# INITIALIZATION
# --------------------------------------------
# Define domain:
L = 2500 # [km]
Dx = 25 # Gridsize [km]
x = np.arange(0, L, Dx)

# Define time array
t = np.arange(0,100)
Dt = 0.1
# Define initial concentration
C = np.zeros((len(t), len(x)))
C[0,46:55] = 1 # Initial concentration

# Define initial wind speed
u0 = 10 # Initial wind speed [m/s]

# --------------------------------------------
# (1) Euler forward in time and upwind in space
# --------------------------------------------
for nt in range(len(t)-1):
    for nx in range(1,len(x)):
        C[nt+1, nx] = -u0*(C[nt,nx]-C[nt,nx-1])/ Dx * Dt+ C[nt,nx]

plt.figure()
plt.plot(x,C[50,:])
plt.show()
# --------------------------------------------
# (2) "Lax-Wendroff" scheme
# --------------------------------------------


# --------------------------------------------
# (3) The spectral method
# --------------------------------------------


#%% animation
    
fig=plt.figure(1)
ax = plt.axes(xlim=(0,L), ylim=(0, 2))
line, = ax.plot([], [], lw=2) #create object to store the data

def init(): #function for initial frame (empty)
    line.set_data([], [])
    return line, 

def anim(i): #function for the time evolution
    line.set_data(x, C[i,:])
    return line, 

#ax.set_title('$\Phi$ for timestep {:.2f} '.format(t)+'$\Lambda$={:.1f} and '.format(Lambda)+'$r$={:.4f}'.format(r))
animat = animation.FuncAnimation(fig,anim, init_func=init, interval=1, blit=True)

