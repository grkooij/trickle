# trickle
Trickle is a 2D Godunov-type fluid simulation that employs the Roe Riemann solver to solver the Euler equations. It uses flat reconstruction and adaptive time-stepping.

## How to use
By default the initial conditions create a Kelvin-Helmholtz instability setup. Define your own initial conditions in `initialise.py`. Make sure pressure is always positive and density is never negative.  After this, define your simulation properties in `constants.py`. Run the simulation with `python trickle.py`
