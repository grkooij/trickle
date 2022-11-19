import numba as nb
from numba import njit
import numpy as np
from constants import RHO, VX, VY, P
from constants import N_COM, GAMMA
from constants import N_X, N_Y, L_X, L_Y
from constants import CFL, CFL_VAR


@njit
def riemann(data, dt):

    directions = [0,1]
    n = [N_X,N_Y]
    
    for direction in directions:

        rh = np.zeros((N_COM, N_X, N_Y), dtype=np.float64)
        
        for i_x in nb.prange(n[direction]):

            p = np.zeros(n[direction-1], dtype=np.float64)
            fl = np.zeros((N_COM, n[direction-1]), dtype=np.float64)  

            smax = 0.
            delta = 1.e-6
            eta = np.zeros(N_COM, dtype=np.float64)
            lam = np.zeros(N_COM, dtype=np.float64)
            rc = np.zeros((N_COM, N_COM), dtype=np.float64)
            v_right = np.zeros((N_COM, n[direction-1]), dtype=np.float64)
            ds = [L_X/N_X, L_Y/N_Y]

            v_n, v_t, v_left = cycle_index(i_x, data, direction)
   
            for v in range(N_COM):
                v_right[v,:] = np.roll(v_left[v,:], -1)

            u_left = prim_to_cons(v_left)
            u_right = prim_to_cons(v_right)

            cs_left = cs(v_left[RHO], v_left[P])
            cs_right = cs(v_right[RHO], v_right[P])

            flux_left = hd_flux(v_left, u_left, n[direction-1], v_n)
            flux_right = hd_flux(v_right, u_right, n[direction-1], v_n)

            for j_x in nb.prange(n[direction-1]):

                u_avg = np.zeros((N_COM), dtype=np.float64)

                vr = v_right[:,j_x]
                vl = v_left[:,j_x]

                dv = vr - vl

                # Roe
                c1 = np.sqrt(vr[RHO]/vl[RHO])
                c2 = 1./(1. + c1)
                c3 = 1. - c2

                u_avg[RHO] = c1*vl[RHO]
                u_avg[VX]  = c2*vl[VX] + c3*vr[VX]
                u_avg[VY]  = c2*vl[VY] + c3*vr[VY]

                v2 = u_avg[VX]*u_avg[VX] + u_avg[VY]*u_avg[VY]

                h_left = 0.5*(vl[VX]*vl[VX] + vl[VY]*vl[VY]) + cs_left[j_x]/(GAMMA-1.)
                h_right = 0.5*(vr[VX]*vr[VX] + vr[VY]*vr[VY]) + cs_right[j_x]/(GAMMA-1.)
                h = c2*h_left + c3*h_right
    
                a = np.sqrt((GAMMA-1.)*(h - 0.5*v2))
  
                smax = max(smax, u_avg[v_n] + a)

                i = 0
                lam[i] = u_avg[v_n] - a
                eta[i] = .5/a/a*(dv[P] - dv[v_n]*u_avg[RHO]*a)
                rc[RHO][i] = 1.
                rc[v_n][i] = u_avg[v_n] - a
                rc[v_t][i] = u_avg[v_t]
                rc[P][i] = h - u_avg[v_n]*a

                i = 1
                lam[i] = u_avg[v_n] + a      
                eta[i] = 0.5/a/a*(dv[P] + dv[v_n]*u_avg[RHO]*a)
                rc[RHO][i] = 1.0
                rc[v_n][i] = u_avg[v_n] + a
                rc[v_t][i] = u_avg[v_t]  
                rc[P][i] = h + u_avg[v_n]*a

                i = 2
                lam[i] = u_avg[v_n]
                eta[i] = dv[RHO] - dv[P]/a/a
                rc[RHO][i] = 1.0
                rc[VX][i] = u_avg[VX]
                rc[VY][i] = u_avg[VY]
                rc[P][i]= 0.5*v2

                i = 3
                lam[i] = u_avg[v_n]
                eta[i]    = u_avg[RHO]*dv[v_t]
                rc[v_t][i] = 1.0
                rc[P][i] = u_avg[v_t]; 

                alam = np.abs(lam)

                # Roe entropy fix
                if alam[0] <= delta:
                    alam[0] = 0.5*lam[0]*lam[0]/delta + 0.5*delta
                if alam[1] <= delta:
                    alam[1] = 0.5*lam[1]*lam[1]/delta + 0.5*delta

                
                flux = flux_left[:,j_x] + flux_right[:,j_x]

                for k in range(N_COM):
                    for l in range(N_COM):
                        flux[k] -= alam[l]*eta[l]*rc[k,l]
                flux *= 0.5    
                p[j_x] = .5*(vl[P] + vr[P])

                for k in range(N_COM):
                    fl[k, j_x] = flux[k]
            
            if direction==0:
                rh[:, i_x, :] = rhs(p, fl, dt, ds[direction-1], n, direction)
            else:
                rh[:, :, i_x] = rhs(p, fl, dt, ds[direction-1], n, direction)

        dt = next_dt(dt, ds[direction-1], smax)

        data = integrate(data, rh)
        if np.isnan(data).any():
            print("nans detected")
    
    return data, dt

@njit
def cs(rho, p):
    gamma = 5./3.

    cs = gamma*p/rho
    return cs

@njit
def hd_flux(v, u, n, v_n):
    flux = np.zeros((N_COM, n), dtype=np.float64)

    flux[RHO] = u[v_n]
    flux[VX]  = u[VX]*v[v_n]
    flux[VY]  = u[VY]*v[v_n]
    flux[P]   = (u[P] + v[P])*v[v_n]

    return flux

@njit
def next_dt(dt, ds, smax):
    

    tmax = CFL*ds/smax
    tnew = CFL_VAR*dt

    next_dt = tnew if dt < tmax and tnew < tmax else min(dt, tmax)
    return next_dt


@njit
def rhs(p, flux, dt, ds, n, direction):
    rhs = np.zeros((N_COM, len(p)), dtype=np.float64)
    v_n = 1 if direction==0 else 2

    flux[v_n,:] += p

    d = dt/ds

    for i in range(n[direction-1]):
        rhs[:,i] = -d*(flux[:, i] - flux[:, i-1])

    return rhs

@njit
def cycle_index(i_x, data, direction):
    v_n = 1 if direction==0 else 2
    v_t = 2 if direction==0 else 1
    if direction==0:
        v_left = data[:, i_x, :]
    else:
        v_left = data[:, :, i_x]

    return v_n, v_t, v_left

@njit
def prim_to_cons(prim):

    cons = np.copy(prim)
    rho = prim[RHO]
    cons[RHO] = rho
    cons[VX] = rho*prim[VX]
    cons[VY] = rho*prim[VY]
    v2 = prim[VX]*prim[VX] + prim[VY]*prim[VY]
    cons[P] = .5*rho*v2 + prim[P]/(GAMMA - 1.)

    return cons

@njit
def cons_to_prim(cons):

    prim = np.copy(cons)

    p2 = cons[VX]*cons[VX] + cons[VY]*cons[VY]
    prim[RHO] = cons[RHO]
    c1 = 1./cons[RHO]
    prim[VX] = c1*cons[VX]
    prim[VY] = c1*cons[VY]
    kin = .5*c1*p2
    prim[P] = (GAMMA-1.)*(cons[P] - kin)

    return prim

@njit
def integrate(data, rhs):

    w0 = .5
    wc = .5

    uc = prim_to_cons(data)
    u0 = np.copy(uc)

    uc += rhs
    uc = w0*u0 + wc*uc
    return cons_to_prim(uc)