

# run simulation to get particle positions
include("DoldMono.jl")
using LaTeXStrings
S = 0.01;
x_o,y_o,t,x,w,xdev,ω,g,k = monoSim(S);

# what does classical theory say

xTheory = copy(x_o);
yTheory = copy(y_o);
for i ∈ 1:length(t), j ∈ 1:length(x)
xTheory[i,j] = w[j] .- (S/k).*cos.(k.*w[j].-sqrt(g*k).*t[i]) .+ 
          sqrt(g/k).* S^2 .* t[i];
yTheory[i,j] = (S/k).*sin.(k.*w[j] .- sqrt(g*k).*t[i]) .+
            0.5*S^2 / k ;
end

# Mean square error
MSE_x = sum([(x_o[i,j]-xTheory[i,j])^2 for i=1:length(t), j=1:length(x)]) / (length(x)*length(t));
MSE_y = sum([(y_o[i,j]-yTheory[i,j])^2 for i=1:length(t), j=1:length(x)]) / (length(x)*length(t));

@gif for i ∈ eachindex(t)
    time = @sprintf("%0.1f",t[i]);
    plot(x_o[i,:],y_o[i,:],xlabel="x (m)",
    ylabel="η (m)",title="S = $S, time = $time (s)",
    xlims=(minimum(x_o),maximum(x_o)),
    ylims = (minimum(y_o),maximum(y_o)),
    label="Numerics")
    plot!(xTheory[i,:],yTheory[i,:],
    label = "Classical Prediction")
end

function MSEtheory(S_range)
    MSE_x = [];
    MSE_y = [];

    for S ∈ S_range
        x_o,y_o,t,x,w,xdev,ω,g,k = monoSim(S);

        # what does classical theory say

        xTheory = copy(x_o);
        yTheory = copy(y_o);
        for i ∈ 1:length(t), j ∈ 1:length(x)
        xTheory[i,j] = w[j] .- (S/k).*cos.(k.*w[j].-sqrt(g*k).*t[i]) .+ 
                sqrt(g/k).* S^2 .* t[i];
        yTheory[i,j] = (S/k).*sin.(k.*w[j] .- sqrt(g*k).*t[i]) .+
                    0.5*S^2 / k ;
        end

        # Mean square error
        sqErrx = sum([(x_o[i,j]-xTheory[i,j])^2 for i=1:length(t), j=1:length(x)]) / (length(x)*length(t));
        sqErry = sum([(y_o[i,j]-yTheory[i,j])^2 for i=1:length(t), j=1:length(x)]) / (length(x)*length(t));
        push!(MSE_x,sqErrx);
        push!(MSE_y,sqErry);
    end
    return MSE_x,MSE_y
end

S_range = range(0.01,0.28,20);
MSE_x,MSE_y = MSEtheory(S_range);

plot(S_range,MSE_x,xlabel=L"S",ylabel = L"\mathrm{MSE} \quad(m^2)",
title="Mean Square Error Averaged over 1 Period",
label = "MSE in x",framestyle= :box)
plot!(S_range,MSE_y,label="MSE in y")
savefig("../../Plots/MSE")