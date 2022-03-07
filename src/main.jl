include("dynProg.jl")
include("plotting.jl")

function instantaneousCost(u, s, p)
    if (u <= s)
        return p*u
    else
        return -Inf
    end
end

function finalCost_aux(s, Pf)
    return Pf*s
end

function dynamics_aux(s, u, w, Smax)
    return round(Int, max(min(s - u + w, Smax), 0))
end

function dependance_noise!(args::Arguments, dhbellman::BellmanFunctions{DH},
    hdbellman::BellmanFunctions{HD}, step::Int64, lastMaxNoise::Int64, multiplier=1)
    y1 = FloatInt[]
    y2 = FloatInt[]
    
    Smax = args.maxStock
    x_axis = last(args.noise.noiseRange):step:lastMaxNoise
    
    for w in x_axis
        println(w)
        
        newArgs = Arguments(args, w, multiplier)

        @time fillvalues!(newArgs, dhbellman)
        @time fillvalues!(newArgs, hdbellman)
        
        hd_middle_value = hdbellman[round(Int, Smax/2),1]
        dh_middle_value = dhbellman[round(Int, Smax/2),1]
    
        push!(y1, dh_middle_value)
        push!(y2, hd_middle_value)
    end
    
    return x_axis, y1, y2
end

function mainplot(args, dhbellman, hdbellman)
    policies_dh = fillvalues!(args, dhbellman)
    policies_hd = fillvalues!(args, hdbellman)

    plotting(args, dhbellman, hdbellman)
    #println(policies_dh)
    #println(policies_hd)
end

function maindepnoise(args, dhbellman, hdbellman, multiplier)
    step = 20
    lastmaxnoise = 300
    xaxis, y1, y2 = dependance_noise!(args, dhbellman, hdbellman, step, lastmaxnoise, multiplier)

    yaxis = [y1, y2]
    
    gr()
    
    display( plot(
    xaxis, yaxis, label=["case 1" "case 2"], legend = :bottomleft, legendfontsize=10,
    line=[:auto :auto],
    color=["blue" "red"], alpha=[0.9 0.9], xlabel="Max Noise "*L"\overline{W}"*" "*L"T="*"$(args.horizon)", xlabelfontsize=16,
    ylabel="Relative distance", ylabelfontsize=16,
    legendtitle="Cases", legendtitlefontsize=12,
    linewidth=3, thickness_scaling = 1, framestyle = :origin
    )
    )
    
    savefig("img/dependencenoise-$(name(args))-$(step)-$(lastmaxnoise).png")
end

function main()
    maxStock = 100
    maxControl = 80
    minControl = 0
    stepControl = 1
    horizon = 8

    minNoise = 0
    maxNoise = 80
    multiplier = 0.001
    noise = Noise{QuadraticDistribution}(minNoise, maxNoise, multiplier)

    p = 2
    P = 10

    dynamics(s, u, w) = dynamics_aux(s, u, w, maxStock)
    instCost(u, s) = instantaneousCost(u, s, p)
    finalCost(s) = finalCost_aux(s, P)

    args = Arguments(maxStock, maxControl, minControl, stepControl, noise, dynamics, instCost, finalCost, horizon)

    dhbellman = BellmanFunctions{DH}(args)
    hdbellman = BellmanFunctions{HD}(args)

    #mainplot(args, dhbellman, hdbellman)
    maindepnoise(args, dhbellman, hdbellman, multiplier)
    
end

main()