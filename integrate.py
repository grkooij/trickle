from numba import njit
import numpy as np
from riemann import riemann
from riemann import prim_to_cons, cons_to_prim

def evolve(data, dt):
    data = riemann(data, dt)
            
    return data



