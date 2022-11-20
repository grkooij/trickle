# trickle
Trickle is a 2D Godunov-type fluid simulation that employs the Roe Riemann solver to solve the Euler equations. It uses flat reconstruction and adaptive time-stepping. Boundary conditions are periodic by default. 

## How to use
By default the initial conditions create a Kelvin-Helmholtz instability setup. Define your own initial conditions in `initialise.py`. Make sure pressure is always positive and density is never negative.  After this, define your simulation properties in `constants.py`. Run the simulation with `python trickle.py`. If you experience issues with NaN's in the simulation domain, it's likely that the CFL variable should be lowered.

## Example
Example snapshot result for the 2D KH instability:
![kh](https://github.com/grkooij/trickle/blob/main/examples/kh.png)