using Plots
using LaTeXStrings
using Printf
include("structs.jl")

function normalize(vec, N)
    if N != 0
        return vec/N
    end
    return vec
end;

function plotting(
    args::Arguments,
    dh_bellman::BellmanFunctions{DH},
    hd_bellman::BellmanFunctions{HD},
    )
    T = args.horizon
    Smax = args.maxStock
    
    for t in T:-1:0
        x_axis = 0:Smax
        y_axis = [
        normalize(view(dh_bellman, :, t+1), (T-t)),
        normalize(view(hd_bellman,:, t+1), (T-t)),
        ]
        
        display( plot(
        x_axis, y_axis, label=["DH" "HD"], legend = :bottomright, legendfontsize=14,
        line=[:auto :auto], linewidht=[20 4],
        color=["blue" "red"], alpha=[0.9 0.9], xlabel="Stock "*L"(s)", xlabelfontsize=16,
        ylabel="Normalized cost " * latexstring("V_{$t}(s)"), ylabelfontsize=16,
        legendtitle=latexstring("t = {$t}"), legendtitlefontsize=16,
        linewidth=3, thickness_scaling = 1, framestyle = :origin
        #title="Normalized cost at time t = $t", titlefontsize=18
        )
        )
        
        #savefig("value_functions_T-$(T-t).png")
    end
end