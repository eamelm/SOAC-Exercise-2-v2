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
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# --------------------------------------------
# INITIALIZATION
# --------------------------------------------
# Define domain:
L = 2500 # [km]
Dx = 25 # Gridsize [km]
x = np.arange(0, L, Dx)

# Define time array
t = np.arange(0,100)

# Define initial concentration
C = np.zeros((len(t), len(x)))
C[46:55] = 1 # Initial concentration

# Define initial wind speed
u0 = 10 # Initial wind speed [m/s]

# --------------------------------------------
# (1) Euler forward in time and upwind in space
# --------------------------------------------
for nt in range(len(t)):
    for nx in range(len(x)):
        C[nx, nt+1] - 

# --------------------------------------------
# (2) "Lax-Wendroff" scheme
# --------------------------------------------


# --------------------------------------------
# (3) The spectral method
# --------------------------------------------

