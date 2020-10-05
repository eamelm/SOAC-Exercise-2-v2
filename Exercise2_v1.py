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

# Define concentration
C0 = 1 # Initial concentration
C = np.zeros(len(x))
C[46:55] = C0

# --------------------------------------------
# (1) Euler forward in time and upwind in space
# --------------------------------------------


# --------------------------------------------
# (2) "Lax-Wendroff" scheme
# --------------------------------------------


# --------------------------------------------
# (3) The spectral method
# --------------------------------------------

