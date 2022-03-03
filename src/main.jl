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
                                hdbellman::BellmanFunctions{HD}, step::Int64)
    y1 = Real[]
    y2 = Real[]
    
    Smax = args.maxStock
    x_axis = args.noise.minNoise:step:args.noise.maxNoise
    
    for w in x_axis
        println(w)
        
        newArgs = Arguments(args, w)

        @time fillvalues!(newArgs, dhbellman)
        @time fillvalues!(newArgs, hdbellman)
        
        hd_middle_value = hdbellman[round(Int, Smax/2),1]
        dh_middle_value = dhbellman[round(Int, Smax/2),1]
    
        push!(y1, dh_middle_value)
        push!(y2, hd_middle_value)
    end
    
    return x_axis, y1, y2
end



function main()
    maxStock = 100
    maxControl = 80
    minControl = 0
    stepControl = 1
    horizon = 8

    minNoise = 0
    maxNoise = 80
    multiplier = 10000
    noise = Noise{LinearDistribution}(minNoise, maxNoise, multiplier)

    p = 2
    P = 10

    dynamics(s, u, w) = dynamics_aux(s, u, w, maxStock)
    instCost(u, s) = instantaneousCost(u, s, p)
    finalCost(s) = finalCost_aux(s, P)

    args = Arguments(maxStock, maxControl, minControl, stepControl, noise, dynamics, instCost, finalCost, horizon)

    dhbellman = BellmanFunctions{DH}(args)
    hdbellman = BellmanFunctions{HD}(args)

    policies_dh = fillvalues!(args, dhbellman)
    policies_hd = fillvalues!(args, hdbellman)

    plotting(args, dhbellman, hdbellman)
    #println(policies_dh)
    #println(policies_hd)
    
end

main()