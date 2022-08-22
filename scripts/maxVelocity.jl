#= 
This script estimates the x-velocity of the particles
and tries to see at which y-positions they occupy.

This is done to try and find where on the wave does most
of the drift occur.

Written by: Aidan Blaser (ablaser@ucsd.edu)
Last Edited: 8/22/2022
=#

# Make sure you run DoldSim.jl first to get simulation
# want x_o and y_o

u = diff(x_o,di)./(t[2]-t[1]);