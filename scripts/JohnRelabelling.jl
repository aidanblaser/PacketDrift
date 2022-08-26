using DrWatson
@quickactivate "PacketDrift"
using ForwardDiff

# Run simulation (adjust parameters there)
include("DoldSim.jl")

# By method of John, want to solve
# z̈ + ig = i r z_w
# first get complex positions of particles
z = (x_o .+ im.*y_o).*L;

# estimate ż and z̈
ż = zeros(length(t)-1,length(z[1,:])) .+ im;
for i ∈ 1:(length(t)-1)
    ż[i,:] = (z[i+1,:]-z[i,:])/(t[i+1]-t[i])./T;
end

z̈ = zeros(length(t)-2,length(z[1,:])) .+ im;
for i ∈ 1:(length(t)-2)
    z̈[i,:] = (ż[i+1,:]-ż[i,:])/(t[i+1]-t[i])./T;
end

# estimate z_w? Assume w = x(t=0)
zw = copy(z);
for i ∈ 1:length(x)-1
    zw[:,i] = (z[:,i+1]) 


plot(z[200,:])
quiver!(real(z[200,:]),imag(z[200,:]),quiver=(real(ż[200,:]),imag(ż[200,:])))


@gif for i ∈ 1:1:length(t)-1
    time = round(T.*t[i],digits=3);
    plot(z[i,:],xlims=(49.8,50.25),ylims=(-0.01,0.01),xlabel="x (m)",
        ylabel="η (m)",title = "S = $S , Δ = $D, t = $time s")
    quiver!(real(z[i,:]),imag(z[i,:]),quiver=(real(ż[i,:])./10,imag(ż[i,:])./10))
end 

plot(z̈[200,:])
plot!(ż[1,:])