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
from scipy.fft import fft, ifft
from matplotlib import animation

# --------------------------------------------
# INITIALIZATION
# --------------------------------------------

# Define domain:
L = 2500 # [km]
Dx = 25 # Gridsize [km]
x = np.arange(0, L, Dx)

# Define time array
Dt = 0.1
t = np.arange(0,100,Dt)

# Define initial concentration
C_e = np.zeros((len(t), len(x)))
C_lw = np.zeros((len(t), len(x)))
C_s= np.zeros((len(t), len(x)))

C_e[0,46:55] = 1   # Initial concentration for Euler forward scheme
C_lw[0,46:55] = 1     # Initial concentration for Lax-Wendroff scheme
C_s[0,46:55] = 1     # Initial concentration for Spectral method
# Define initial wind speed
u0 = 10 # Initial wind speed [m/s]

# --------------------------------------------
# (1) Euler forward in time and upwind in space
# --------------------------------------------
for nt in range(len(t)-1):  # loop over time
    for nx in range(1, len(x)):      # loop over space
        C_e[nt+1, nx] = -u0 * (C_e[nt, nx]-C_e[nt, nx-1]) / Dx * Dt + C_e[nt, nx]


# --------------------------------------------
# (2) "Lax-Wendroff" scheme
# --------------------------------------------
for nt in range(len(t)-1):
    for nx in range(1,len(x)-1):
        C_lw[nt+1, nx] = (-u0*(C_lw[nt,nx+1]-C_lw[nt,nx-1])/ (2*Dx) + 
        (u0**2 * Dt/2*(C_lw[nt,nx+1]-2*C_lw[nt,nx]+C_lw[nt,nx-1])/Dx**2)) * Dt+ C_lw[nt,nx]

# --------------------------------------------
# (3) The spectral method
# --------------------------------------------
# N = len(x)
# alpha = np.zeros(int(N/2))
# beta = np.zeros(int(N/2))
# C_k = np.zeros((len(t),int(N/2)),dtype=complex)

# alpha[0] = 1/N * np.sum(C_s[:-1])

# for nt in range(len(t)-1):
#     for k in range(int(N/2)):
#         for j in range(N-1):
#             alpha[k] +=  C_s[nt,j]*np.cos(-2*np.pi*k*j/N)
#             beta[k] += C_s[nt,j] * np.sin(-2*np.pi*k*j/N)
#         alpha[k] *= 2/N
#         C_k[nt,k] = alpha[k] + beta[k] * 1j
#         C_k[nt+1] = -2*np.pi* 1j * u0 * k/L * C_k[nt] * Dt

# C_k_inv = np.zeros((int(N/2),N),dtype=complex)
# for nt in range(len(t)-1):
#     for j in range(N):
#         for k in range(int(N/2)):
#             C_k_inv[k,j] = C_k[nt,k]*np.exp(1j*k*j)
#         C_s[nt,j] = (np.sum(C_k_inv[:,j])/(N/2))
#%%        
C_ft = np.zeros((np.shape(C_s)),dtype=complex)
######### OPTION 2:
# for nt in range(len(t)-1):
C_ft[0,:] = fft(C_s[0,:])
    # C_ft[nt+1,:] = C_ft[nt,:]*Dt
plt.figure()
plt.plot(x,C_ft[0,:])
plt.show() 
#%% animation
    
fig=plt.figure(1)
ax = plt.axes(xlim=(0,L), ylim=(0, 2))
line, = ax.plot([], [], lw=2)
line_lw, = ax.plot([], [], lw=2) #create object to store the data

def init(): #function for initial frame (empty)
    # line.set_data([], [])
    line_lw.set_data([], [])
    return line, line_lw,

def anim(i): #function for the time evolution
    # line.set_data(x, C_e[i,:])
    line_lw.set_data(x,C_s[i,:])
    return line, line_lw,

#ax.set_title('$\Phi$ for timestep {:.2f} '.format(t)+'$\Lambda$={:.1f} and '.format(Lambda)+'$r$={:.4f}'.format(r))
animat = animation.FuncAnimation(fig,anim, init_func=init, interval=1, blit=True)

