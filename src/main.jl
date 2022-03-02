include("dynProg.jl")
include("plotting.jl")

function linearLaw(minNoise::Int64, maxNoise::Int64, multiplier::Real, w)
    minProb = 2/((multiplier+1)*(maxNoise - minNoise + 1))
    deltaProb = (multiplier - 1) * minProb / (maxNoise - minNoise)
    return minProb + deltaProb * (w - minNoise)
end

function uniformLaw(minNoise::Int64, maxNoise::Int64, multiplier::Real, w)
    return 1/(maxNoise - minNoise + 1)
end

function quadraticLaw(minNoise::Int64, maxNoise::Int64, multiplier::Real, w)
    return 6*(w - minNoise)*(maxNoise - w)/(maxNoise - minNoise)^3
end

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

function main()
    maxStock = 100
    maxControl = 80
    minControl = 0
    stepControl = 1
    horizon = 25

    minNoise = 0
    maxNoise = 120
    multiplier = 10000
    law(w) = quadraticLaw(minNoise, maxNoise, multiplier, w)
    noise = Noise(maxNoise, minNoise, law)

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
end

main()