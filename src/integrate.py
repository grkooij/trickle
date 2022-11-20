from numba import njit
import numpy as np
from src.riemann import riemann

def evolve(data, dt):
    data = riemann(data, dt)
            
    return data



