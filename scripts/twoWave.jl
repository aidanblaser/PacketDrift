include("DoldDuo.jl")

S1 = 0.001;
S2 = 0.001;
k1 = 1;
k2 = 4;

x_o,y_o,u,v,t,x,w,xdev,ω1,ω2,g = duoSim(S1,S2,k1=k1,k2=k2);

cd(projectdir())
pyplot()
anim = @animate for i ∈ 1:length(t)
    time = @sprintf("%0.1f",t[i]);
    plot(x_o[i,:] , y_o[i,:],  xlabel="x (m)",ylabel="η (m)",
    title = "S = ($S1,$S2),   k = ($k1,$k2),   t = $time s",label="Numerics",
    xlims=(minimum(x_o),maximum(x_o)),
    ylims=(minimum(y_o),maximum(y_o)),
    framestyle= :box)
    #plot!(x_o[1:i,250].*L,y_o[1:i,250].*L,label="",color="red")
    #scatter!([x_o[i,250].*L],[y_o[i,250].*L],label="")
end 
gif(anim, "Plots/twoWaves.gif", fps = 6)

# Idea: Want to find η(x) so I have to interpolate
# at each time

η = [];
for j ∈ 1:length(t)
push!(η,interpolate((x_o[j,:],),y_o[j,:],Gridded(Linear())));
end

# for a given point in x, find a time series
ind = Int(ceil(length(x)/2));

time_series = Float64[];
for j ∈ 1:length(t)-1
    interp = η[j]
    push!(time_series,interp(x[ind]))
end

plot(t[1:end-1],time_series)

# Take spectrum to find frequencies
Ft = 1/(t[2]-t[1]);

out,freq = positiveFFT(time_series,Ft);

plotlyjs()
plot(freq*2π,abs.(out),xlabel="ω (rad/s)",
ylabel="abs(FFT(y))",
title="Time Spectrum for Two Waves",
framestyle= :box)
