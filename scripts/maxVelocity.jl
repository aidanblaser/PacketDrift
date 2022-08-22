#= 
This script estimates the x-velocity of the particles
and tries to see at which y-positions they occupy.

This is done to try and find where on the wave does most
of the drift occur.

Written by: Aidan Blaser (ablaser@ucsd.edu)
Last Edited: 8/22/2022
=#
using DrWatson
@quickactivate "PacketDrift"
using Plots
using Distributions


# Make sure you run DoldSim.jl first to get simulation
# want x_o and y_o

u = copy(x_o);

for i ∈ 1:(length(t)-1)
    u[i,:] = (x_o[i+1,:]-x_o[i,:])/(t[i+1]-t[i]);
end

@gif for i ∈ 1:2:length(t)
    time = round(T.*t[i],digits=3);
    plot(x_o[i,:].*L, u[i,:].*L./T, xlims=(0,wl),
    ylims = (-0.15,0.15), xlabel="x (m)",ylabel="u (m/s)",
    title = "S = $S , Δ = $D, t = $time s")
end 


histogram2d(u[:],y_o[:])