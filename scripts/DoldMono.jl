# General Dold Simulation for Monochromatic
# Differs from DoldSim as it mostly uses single waves

using DrWatson
@quickactivate "PacketDrift"
using FFTW
using Plots
using DelimitedFiles
using Interpolations
using ForwardDiff

function monoSim(S,k=1)
# Adjustable Parameters
g = 9.81; # acc due to gravity. make sure Dold is set to same value
ω = sqrt(g*k);
tl = 2π / ω ; # length of simulation (seconds)
N = 512; # Number of surface points
ML = 2π*k; # Physical length of channel
x =collect( 0 : ML/N : (ML) * (1 - 1/N)); # domain
BW = 0; # bandwidth (set to 0 for no focusing packet)
wl = ML; # wl parameter is used in Dold code

# Define surface displacement eta. Linear dispersive focusing 
eta1 = zeros(length(x),1);
eta1 = (S/k).*sin.(k.*x);
#filter = exp.(-(x.-5).^2 ./ (5^2));
eta2 = eta1;


## we now want to find the velocity potential
Fs = 1 / ( abs(x[2] - x[1]) ); # define sampling frequency in space
function positiveFFT(x,Fs)
    N=length(x); 
    k=collect(0:N-1); 
    T=N/Fs; 
    freq=k/T; #create the frequency range 
    X=FFTW.fft(x)/N; # normalize the data
    cutOff = Int(ceil(N/2)); 
    X = X[1:cutOff]; 
    freq = freq[1:cutOff];
    return X, freq
end
out,freq = positiveFFT(eta2, Fs); #perform fft
# check fft on eta
eta = zeros(1, length(x)); 
for j ∈ 1 : length(x)
    eta[1,j] = sum( real( 2 .* out[2:end].*
    exp.( im .* ( 2 * π .* freq[2:end].*(x[j]) ) )));
end
# check S for the windowed wave form
eps2 = sum((abs.(2 .*out[2:end]) .* 2 .* π .* freq[2:end] )); #
# now find phi
phi = zeros(1, length(x)); 
for j ∈ 1 : length(x)
    phi[1,j] = sum(imag.(2 .* sqrt(g) .* out[2:end] ./ 
        sqrt.(2 .* π .* freq[2:end]).*
        exp.( im .* (( 2 .* π .* freq[2:end] ).*
        (x[j]) ))));
end
# BW=powerbw(eta,2*pi*Fs,[],10)/k ; # check bandwidth

x_f = x; # shift x-axis so it goes from 0 to wl
y_f = eta;
f_f = phi;

plot( x_f, [y_f',f_f'], xlabel='x', label=["η" "ϕ"],
        legendfontsize=18)
#set(gca, 'fontsize', 22)
#l1 = legend('$\eta$' , '$\phi$');
#set(l1, 'interpreter' , 'latex')


## 
#tl = (xb+20)*D; # length of simulation

# Make sure you're in right directory
cd(projectdir()) 
cd("scripts/dold_files")

rm("xc.txt",force=true)
rm("yc.txt",force=true)
rm("fc.txt",force=true)
rm("bw.txt",force=true)
rm("wl.txt",force=true)
rm("S.txt",force=true)
rm("bw0.txt",force=true)
rm("S0.txt",force=true)
rm("C.txt",force=true)
# delete k.txt
writedlm("xc.txt",x_f)
writedlm("yc.txt", y_f)
writedlm("fc.txt", f_f)
writedlm("S.txt", eps2)
writedlm("S0.txt", S)
#writedlm("bw.txt", BW)
#writedlm("bw0.txt", BW)
#writedlm("C.txt", BW)
writedlm("wl.txt", wl)
writedlm("tl.txt", tl)
# save k.txt k -ascii
run(`./run2.sh`);

# Open files and save as variables
x_o = readdlm("x.txt");
y_o = readdlm("y.txt");
p_o = readdlm("phi.txt");
to = readdlm("t.txt");
t = to[:,1];


# re-dimensionalize variables
L = wl / 2 / π; #length scale
T = sqrt( wl / 2 / π ); # time scale
xi = 0 : wl / N : wl - wl / N; # spatial grid

x_o *= L;
y_o *= L;
t = t*T;

# Estimate velocity
u = copy(x_o);
v = copy(y_o);
for i ∈ 1:length(x), j ∈ 1:length(t)
    xint = interpolate((t,),x_o[:,i],Gridded(Linear()));
    u[j,i] = only(gradient(xint,t[j]));
    yint = interpolate((t,),y_o[:,i],Gridded(Linear()));
    v[j,i] = only(gradient(yint,t[j]));
end

#= look at the data
pyplot()
anim = @animate for i ∈ 1:length(t)
    time = round(t[i],digits=3);
    plot(x_o[i,:] , y_o[i,:],  xlabel="x (m)",ylabel="η (m)",
    title = "S = $S, t = $time s")
    quiver!(x_o[i,1:10:end],y_o[i,1:10:end],
    quiver=(u[i,1:10:end],v[i,1:10:end]))
    #plot!(x_o[1:i,250].*L,y_o[1:i,250].*L,label="",color="red")
    #scatter!([x_o[i,250].*L],[y_o[i,250].*L],label="")
end 
gif(anim, "anim_fps15.gif", fps = 6)
=#

w = x_o[1,:];
xdev = copy(x_o);
for i ∈ 1:length(x)
    xdev[:,i] = x_o[:,i] .- w[i]
end

k = ω^2 / g;

return x_o,y_o,t,x,w,xdev,ω,g,k
end