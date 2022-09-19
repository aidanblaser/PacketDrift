#=
Goal of this code is to be able to numerically
compute r(w,t) for the flow
=#

# run simulation to get particle locations
include("DoldMono.jl")

S = 0.01;
x_o,y_o,u,v,t,x,w,xdev,ω,g,k = monoSim(S);

# Estimate acceleration
ax = copy(x_o);
ay = copy(y_o);
zw = copy(x_o).+ im;
for i ∈ 1:length(x), j ∈ 1:length(t)
    uint = interpolate((t,),u[:,i],Gridded(Linear()));
    ax[j,i] = only(gradient(uint,t[j]));
    vint = interpolate((t,),v[:,i],Gridded(Linear()));
    ay[j,i] = only(gradient(vint,t[j]));
    zwint = interpolate((w,),x_o[j,:] .+ im.*y_o[j,:],Gridded(Linear()));
    zw[j,i] = only(gradient(zwint,w[i]))
end

#=
cd(projectdir())
pyplot()
anim = @animate for i ∈ 1:length(t)
    time = @sprintf("%0.1f",t[i]);
    plot(x_o[i,:] , y_o[i,:],  xlabel="x (m)",ylabel="η (m)",
    title = "S = $S, t = $time s",label="Numerics",
    xlims=(minimum(x_o),maximum(x_o)),
    ylims=(minimum(y_o.-1.5.*ay),maximum(y_o.+1.5*ay)),
    framestyle= :box)
    quiver!(x_o[i,1:10:end],y_o[i,1:10:end],
    quiver=(ax[i,1:10:end],ay[i,1:10:end]),
    label="Velocity")
    #plot!(x_o[1:i,250].*L,y_o[1:i,250].*L,label="",color="red")
    #scatter!([x_o[i,250].*L],[y_o[i,250].*L],label="")
end 
gif(anim, "Plots/acceleration.gif", fps = 6)
=#

# Last acceleration time is weird, ignore for now
r = zeros(length(t)-1,length(x)).+ im;
for i ∈ 1:length(t)-1
    r[i,:] = (ax[i,:] .+ im.* ay[i,:]) ./ (im.*zw[i,:]);
end 

plot(r[1,:])
plot(r[2,:])
plot(r[12,:])
plot!(r[18,:])
plot!(r[21,:])