import numpy as np
from constants import N_X, N_Y, L_X, L_Y
import matplotlib.pyplot as plt

def initialise():

    rho = np.zeros((N_X, N_Y), dtype=np.float64)
    vx = np.zeros((N_X, N_Y), dtype=np.float64)
    vy = np.zeros((N_X, N_Y), dtype=np.float64)
    p = np.zeros((N_X, N_Y), dtype=np.float64)

    vel = 1.0

    x = np.linspace(0, L_X, N_X, dtype=np.float64)
    y = np.linspace(0, L_Y, N_Y, dtype=np.float64)

    xx, yy = np.meshgrid(x, y, indexing='ij')
    b_upper = 0.7
    b_lower = 0.3
    amp = 0.01

    for i in range(N_X):
        for j in range(N_Y):

            Y = yy[i,j]
            X = xx[i,j]

            if X > b_lower and X < b_upper:
                vx[i,j]  = -vel
                rho[i,j] = 2.0
            
            if X <= b_lower:
                vx[i,j]  = vel
                rho[i,j] = 1.0

            if X >= b_upper:
                vx[i,j]  = vel
                rho[i,j] = 1.0

            sigma = 0.05/np.sqrt(2)
            vy[i,j] = amp*np.sin(2*np.pi*Y)*(np.exp(-(X - b_lower)*(X - b_lower)/2./sigma/sigma) + np.exp(-(X - b_upper)*(X - b_upper)/2./sigma/sigma))

    p += 2.5

    data = np.array([rho, vx, vy, p])

    return data