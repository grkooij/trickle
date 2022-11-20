from numba import njit
import numpy as np
from riemann import riemann

def evolve(data, dt):
    data = riemann(data, dt)
            
    return data



