

# run simulation to get particle positions
include("DoldMono.jl")
using LaTeXStrings
using Printf
S = 0.01;
x_o,y_o,u,v,t,x,w,xdev,ω,g,k = monoSim(S);

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

anim2 = @animate for i ∈ eachindex(t)
    time = @sprintf("%0.1f",t[i]);
    plot(x_o[i,:],y_o[i,:],xlabel="x (m)",
    ylabel="η (m)",title="S = $S, time = $time (s)",
    xlims=(minimum(x_o),maximum(x_o)),
    ylims = (minimum(y_o),maximum(y_o)),
    label="Numerics",framestyle=:box)
    plot!(xTheory[i,:],yTheory[i,:],
    label = "Classical Prediction")
end
cd(projectdir())
gif(anim2,"Plots/Numerics_v_Theory.gif",fps=6)

function MSEtheory(S_range)
    MSE_x = Vector{Type:Float64};
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
label = "MSE in x",framestyle= :box,
legend = :topleft)
plot!(S_range,MSE_y,label="MSE in y")
savefig("../../Plots/MSE")

#Try a least squares
# S^2 fit 
S2 = S_range.^2;
d = convert(Vector{Float64},MSE_x);
α = S2 \ d;

plot!(S_range,α .*S2,label=L"\alpha S^2")
cd(projectdir())
savefig("Plots/SquareFit")

# S^2 + S^4
S4 = S_range.^4;
SMat = hcat(S2,S4);
β = SMat \ d;
γ = S4 \ d;
beta = round.(β,digits=2);
plot!(S_range,β[1] .*S2 .+ β[2] .*S4,label=
string("$(beta[1]) ",L"S^2 + "," $(beta[2]) ",L"S^4"))
cd(projectdir())
savefig("Plots/QuarFit")

plot!(S_range,γ.*S4)