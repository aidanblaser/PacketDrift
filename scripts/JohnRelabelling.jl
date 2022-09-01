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
zw = zeros(length(t),length(x)-1) .+ im;
for i ∈ 1:length(x)-1
    zw[:,i] = (z[:,i+1]- z[:,i])./ (x[i+1]-x[i])
end 

# estimate r(w,t)
r = (z̈[:,1:end-1] .+ 0 .*im*g)./(im.*zw[1:end-2,:]);

plot(r[2,:])

@gif for i ∈ 1:length(t)-2
    time = round(T.*t[i],digits=3);
    plot(r[i,:],xlims=(9.4-g,10.25-g),ylims=(-0.4,0.4),
    xlabel="Re[r(w,t)]",ylabel="Im[r(w,t)]",
    title="S = $S, Δ = $D, t = $time s",
    label = "r(w,t)")
end


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