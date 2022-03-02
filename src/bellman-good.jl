using Plots
using LaTeXStrings
using Printf

struct Arguments
    maxStock::Int64
    maxControl::Int64
    minControl::Int64
    stepControl::Int64
    maxNoise::Int64
    minNoise::Int64
    instantaneousCostConstant::Number
    finalCostConstant::Number
    horizon::Int64
end

struct Noise
    maxNoise::Int64
    minNoise::Int64
    law::Function
end

struct Variables
    maxStock::Int64
    maxControl::Int64
    minControl::Int64
    stepControl::Int64
    noise::Noise
    instantaneousCost::Function
    finalCost::Function
    horizon::Int64
end


function Arguments(args::Arguments, newMaxNoise::Int64)
    newArgs = Arguments(args.maxStock, args.maxControl, 0, 1, newMaxNoise, 0,
    args.instantaneousCostConstant, args.finalCostConstant, args.horizon,
    )
    return newArgs
end

function instantaneousCost2(u, s, w, args::Arguments)
    Smax = args.maxStock
    p = args.instantaneousCostConstant
    if (u <= min(s+w, Smax))
        return p*u
    else
        return -Inf
    end
end

function instantaneousCost(u, s, w, args::Arguments)
    Smax = args.maxStock
    p = args.instantaneousCostConstant
    if (u <= s)
        return p*u
    else
        return -Inf
    end
end

function finalCost(s, args::Arguments)
    Pf = args.finalCostConstant
    return Pf*s
end

function dynamics_aux(s, u, w, args::Arguments)
    Smax = args.maxStock
    return round(Int, max(min(s - u + w, Smax), 0))
end

function fillvalues!(args::Arguments, bellman_function::Array{U, 2}) where U <: Real
    Smax = args.maxStock
    Umax = args.maxControl
    Umin = args.minControl
    Wmax = args.maxNoise
    Wmin = args.minNoise
    Ustep = args.stepControl
    
    Winv = 1/(Wmax - Wmin)
    
    T = args.horizon
    K(s) = finalCost(s, args)
    L(u,s,w) = instantaneousCost(u, s, w, args)
    dynamics(s, u, w) = dynamics_aux(s, u, w, args)
    
    policies = zeros(T)

    for s in 1:Smax+1
        bellman_function[s,T+1] = K(s - 1)
    end
    
    for t in T-1:-1:0, s in 0:Smax # t = T-1, ..., 0, # s = 0, 1, ..., Smax
        max_expected = -Inf # needs to be -Inf
        for u in Umin:Ustep:Umax # u = Umin, Umin+1, ..., Umax 
            expected_value = 0
            for i in Wmin:Wmax # i = 0, 1, ..., Wmax
                w = Wmin + i * Winv * (Wmax - Wmin)
                realization = L(u, s, w) + bellman_function[dynamics(s, u, w) + 1, t+2]
                expected_value += realization 
            end
            expected_value /= (Wmax - Wmin + 1) # Wmax + 1 realizations
            if max_expected < expected_value
                max_expected = expected_value #max(expected_value, max_expected)
                policies[t+1] = u
            end
        end
        bellman_function[s+1,t+1] = max_expected
    end

    return policies
end

function fillvalues_hd!(args::Arguments, hd_bellman_function::Array{U, 2}) where U <: Real
    Smax = args.maxStock
    Umax = args.maxControl
    Umin = args.minControl
    Wmax = args.maxNoise
    Wmin = args.minNoise
    Ustep = args.stepControl
    
    Winv = 1/(Wmax - Wmin)
    
    T = args.horizon
    K(s) = finalCost(s, args)
    L(u,s,w) = instantaneousCost(u, s, w, args)
    
    dynamics(s, u, w) = dynamics_aux(s, u, w, args)
    
    for s in 1:Smax+1
        hd_bellman_function[s,T+1] = K(s - 1)
    end
    
    for t in T-1: -1: 0 # t = T-1, ..., 0
        for s in 0:Smax # s = 0, 1, ..., Smax
            expected_value = 0
            for i in 0:Wmax # i = 0, 1, ..., Nw
                w = Wmin + i * Winv * (Wmax - Wmin)
                max_value = -Inf
                for u in Umin:Ustep:Umax # u = Umin, ..., Umax 
                    max_value = max(max_value, L(u, s, w) + hd_bellman_function[dynamics(s, u, w)+1, t+2])
                end
                expected_value += max_value
            end
            expected_value /= (Wmax + 1)
            hd_bellman_function[s+1,t+1] = expected_value
        end
    end
end

function fillvalues_hd_2!(args::Arguments, hd_bellman_function::Array{U, 2}) where U <: Real
    Smax = args.maxStock
    Umax = args.maxControl
    Umin = args.minControl
    Wmax = args.maxNoise
    Wmin = args.minNoise
    Ustep = args.stepControl
    
    Winv = 1/(Wmax - Wmin)
    
    T = args.horizon
    K(s) = finalCost(s, args)
    L(u,s,w) = instantaneousCost2(u, s, w, args)
    
    dynamics(s, u, w) = dynamics_aux(s, u, w, args)
    
    for s in 1:Smax+1
        hd_bellman_function[s,T+1] = K(s - 1)
    end
    
    for t in T-1: -1: 0 # t = T-1, ..., 0
        for s in 0:Smax # s = 0, 1, ..., Smax
            expected_value = 0
            for i in 0:Wmax # i = 0, 1, ..., Nw
                w = Wmin + i * Winv * (Wmax - Wmin)
                max_value = -Inf
                for u in Umin:Ustep:Umax # u = Umin, ..., Umax 
                    max_value = max(max_value, L(u, s, w) + hd_bellman_function[dynamics(s, u, w)+1, t+2])
                end
                expected_value += max_value
            end
            expected_value /= (Wmax + 1)
            hd_bellman_function[s+1,t+1] = expected_value
        end
    end
end

function showPolicies()
    argsTup = (100, 80, 0, 1, 10, 0, 2, 10000, 25)
    args = Arguments(argsTup...)

    dh_bellman = zeros(args.maxStock + 1, args.horizon + 1)

    policies = fillvalues!(args, dh_bellman)
    for t in 1:length(policies)
        println("Time $(t-1) : decision $(policies[t]).")
    end
end

function simulate(args::Arguments)
    Smax = args.maxStock
    T = args.horizon
    
    dh_bellman_values = zeros(Smax+1, T+1)
    hd_bellman_values = zeros(Smax+1, T+1)
    hd_bellman_values_2 = zeros(Smax+1, T+1)
    
    @time fillvalues!(args, dh_bellman_values)
    @time fillvalues_hd!(args, hd_bellman_values)
    @time fillvalues_hd_2!(args, hd_bellman_values_2)
    
    plotting(args, dh_bellman_values, hd_bellman_values, hd_bellman_values_2)
end

function simulate(
    args::Arguments,
    dh_bellman::Array{U, 2}, 
    hd_bellman::Array{U, 2}, 
    hd_bellman2::Array{U, 2}
    ) where U <: Real
    
    fillvalues!(args, dh_bellman)
    fillvalues_hd!(args, hd_bellman)
    fillvalues_hd_2!(args, hd_bellman2)
end

function oldsimulate(args::Arguments)
    Smax = args.maxStock
    T = args.horizon
    
    dh_bellman_values = zeros(Smax+1, T+1)
    hd_bellman_values = zeros(Smax+1, T+1)
    
    fillvalues!(args, dh_bellman_values)
    fillvalues_hd!(args, hd_bellman_values)
    
    return (hd_bellman_values[round(Int, Smax/2),1] - dh_bellman_values[round(Int, Smax/2),1])/hd_bellman_values[round(Int, Smax/2),1]
end

function oldsimulate2(args::Arguments)
    Smax = args.maxStock
    T = args.horizon
    
    dh_bellman_values = zeros(Smax+1, T+1)
    hd_bellman_values = zeros(Smax+1, T+1)
    
    fillvalues!(args, dh_bellman_values)
    fillvalues_hd_2!(args, hd_bellman_values)
    
    return (hd_bellman_values[round(Int, Smax/2),1] - dh_bellman_values[round(Int, Smax/2),1])/hd_bellman_values[round(Int, Smax/2),1]
end

function oldmultiple_simulations(wmaxmultiplier = 2)
    diff = 0
    smax, umax, wmax, pmax, Pmax = 0,0,0,0,0
    umin = 0
    x = 0
    wmin = 0
    ustep = 20
    sstart = 50
    horizon = 8
    
    #p = 10
    #P = 5
    
    for s in sstart:50:300, u in 0.5*s:ustep:0.9*s, w in 0.2*s:ustep:4*s, p in 2:4:10, P in 5:5:25
        #w = wmaxmultiplier*s
        println("s, u, w, p, P : $(s), $(u), $(w), $(p), $(P)")
        args = Arguments(s,u,umin,ustep,w, wmin, p, P, horizon)
        x = @time oldsimulate2(args)
        if (diff < x)
            diff = x
            smax, umax, wmax, pmax, Pmax = s, u, w, p, P
            printstyled("Found diff $(x) for s, u, w, p, P = $(s), $(u), $(w), $p, $P\n"; color = :green )
        else
            printstyled("Found diff $(x) for s, u, w, p, P = $(s), $(u), $(w), $p, $P "; color = :red )
            printstyled("still $(diff) for s, u, w, p, P = $(smax), $(umax), $(wmax), $pmax, $Pmax\n"; color = :green )
        end
    end
end

function multiple_simulations_twodiff()
    diff1, diff2 = 0, 0
    smax, umax, wmax, pmax, Pmax = 0,0,0,0,0
    umin = 0
    x = 0
    wmin = 0
    ustep = 20
    sstart = 100
    horizon = 8
    
    p = 2
    P = 10
    
    for s in sstart:50:100, u in 0.2*s:ustep:0.9*s, w in 0.2*s:ustep:3*s#, p in 2:4:10, P in 5:5:25
        Smax = s
        T = horizon
        dh_bellman_values = zeros(Smax+1, T+1)
        hd_bellman_values = zeros(Smax+1, T+1)
        hd_bellman_values_2 = zeros(Smax+1, T+1);
        println("s, u, w, p, P : $(s), $(u), $(w), $(p), $(P)")
        args = Arguments(s,u,umin,ustep,w, wmin, p, P, horizon)
        @time simulate(args, dh_bellman_values, hd_bellman_values, hd_bellman_values_2)
        
        x1 = (hd_bellman_values[round(Int, Smax/2),1] - dh_bellman_values[round(Int, Smax/2),1])/hd_bellman_values[round(Int, Smax/2),1]
        x2 = (hd_bellman_values_2[round(Int, Smax/2),1] - dh_bellman_values[round(Int, Smax/2),1])/hd_bellman_values_2[round(Int, Smax/2),1]
        
        if (diff1 <= x1 && diff2 <= x2)
            diff1, diff2 = x1, x2
            smax, umax, wmax, pmax, Pmax = s, u, w, p, P
            printstyled("Found diff $(x1) $x2 for s, u, w, p, P = $(s), $(u), $(w), $p, $P\n"; color = :green )
        else
            printstyled("Found diff $(x1) $x2 for s, u, w, p, P = $(s), $(u), $(w), $p, $P "; color = :red )
            printstyled("still $(diff1) $(diff2) for s, u, w, p, P = $(smax), $(umax), $(wmax), $pmax, $Pmax\n"; color = :green )
        end
    end
end

function multiple_simulations_wdependance()
    diff1, diff2 = 0, 0
    smax, umax, wmax, pmax, Pmax = 0,0,0,0,0
    umin = 0
    x = 0
    wmin = 0
    ustep = 20
    sstart = 100
    horizon = 8
    
    dh_bellman_values = zeros(Smax+1, T+1)
    hd_bellman_values = zeros(Smax+1, T+1)
    hd_bellman_values_2 = zeros(Smax+1, T+1)
    
    s, u, p, P = 100, 80, 2, 10
    
    args = Arguments(s,u,umin,ustep,w, wmin, p, P, horizon)
    @time simulate(args, dh_bellman_values, hd_bellman_values, hd_bellman_values_2)
    
    x1 = (hd_bellman_values[round(Int, Smax/2),1] - dh_bellman_values[round(Int, Smax/2),1])/hd_bellman_values[round(Int, Smax/2),1]
    x2 = (hd_bellman_values_2[round(Int, Smax/2),1] - dh_bellman_values[round(Int, Smax/2),1])/hd_bellman_values_2[round(Int, Smax/2),1]
end

function plotting(
    args::Arguments,
    dh_bellman::Array{U, 2},
    hd_bellman::Array{U, 2},
    hd_bellman2::Array{U, 2}
    ) where U <: Real
    T = args.horizon
    Smax = args.maxStock
    
    for t in T:-1:0
        x_axis = 0:Smax
        y_axis = [
        normalize(view(dh_bellman, :, t+1), (T-t)),
        normalize(view(hd_bellman,:, t+1), (T-t)),
        normalize(view(hd_bellman2,:, t+1), (T-t))
        ]
        
        display( plot(
        x_axis, y_axis, label=["DH" "HD" "HD2"], legend = :bottomright, legendfontsize=14,
        line=[:auto :auto :auto], linewidht=[20 4 1],
        color=["blue" "red" "orange"], alpha=[0.9 0.9 0.9], xlabel="Stock "*L"(s)", xlabelfontsize=16,
        ylabel="Normalized cost " * latexstring("V_{$t}(s)"), ylabelfontsize=16,
        legendtitle=latexstring("t = {$t}"), legendtitlefontsize=16,
        linewidth=3, thickness_scaling = 1, framestyle = :origin
        #title="Normalized cost at time t = $t", titlefontsize=18
        )
        )
        
        #savefig("value_functions_T-$(T-t).png")
    end
end

function normalize(vec, N)
    if N != 0
        return vec/N
    end
    return vec
end;

function olddependance_noise!(
    args::Arguments,
    dh_bellman::Array{U, 2}, 
    hd_bellman::Array{U, 2}, 
    hd_bellman2::Array{U, 2}
    ) where U <: Real
    
    y1 = Real[]
    y2 = Real[]
    
    Smax = args.maxStock
    x_axis = args.minNoise:args.stepControl:args.maxNoise
    
    for w in x_axis
        newArgs = Arguments(args, w)
        println(w)
        
        @time fillvalues!(newArgs, dh_bellman)
        @time fillvalues_hd!(newArgs, hd_bellman)
        @time fillvalues_hd_2!(newArgs, hd_bellman2)
        
        hd_middle_value = hd_bellman[round(Int, Smax/2),1]
        dh_middle_value = dh_bellman[round(Int, Smax/2),1]
        hd_middle_value2 = hd_bellman2[round(Int, Smax/2),1]
        
        relative_distance = (hd_middle_value - dh_middle_value)/hd_middle_value
        relative_distance2 = (hd_middle_value2 - dh_middle_value)/hd_middle_value2
        
        push!(y1, relative_distance)
        push!(y2, relative_distance2)
    end
    
    return x_axis, y1, y2
end

function dependance_noise!(
    args::Arguments,
    dh_bellman::Array{U, 2}, 
    hd_bellman::Array{U, 2}
    ) where U <: Real
    
    y1 = Real[]
    y2 = Real[]
    
    Smax = args.maxStock
    x_axis = args.minNoise:args.stepControl:args.maxNoise
    
    for w in x_axis
        newArgs = Arguments(args, w)
        println(w)
        
        @time fillvalues!(newArgs, dh_bellman)
        @time fillvalues_hd!(newArgs, hd_bellman)
        
        hd_middle_value = hd_bellman[round(Int, Smax/2),1]
        dh_middle_value = dh_bellman[round(Int, Smax/2),1]
    
        push!(y1, dh_middle_value)
        push!(y2, hd_middle_value)
    end
    
    return x_axis, y1, y2
end

function main2()
    maxStock = 100
    maxControl = 80
    minControl = 0
    stepControl = 1
    maxNoise = 180
    minNoise = 0
    instantaneousCostConstant = 10
    finalCostConstant = 5
    horizon = 10
    
    args = Arguments(maxStock, maxControl, minControl, stepControl, maxNoise, minNoise,
    instantaneousCostConstant, finalCostConstant, horizon)
    
    println(oldsimulate(args))
end

function olddependance_noise_test(horizon::Int64)
    noisetest1 = (100, 80, 0, 20, 500, 20, 2, 0, horizon)
    noisetest2 = (200, 160, 0, 20, 500, 20, 6, 25, 4)
    noisetest3 = (50, 45, 0, 20, 500, 20, 10, 5, 4)
    noisetest4 = (150, 130, 0, 20, 500, 20, 2, 10, 4)
    noisetest5 = (150, 130, 0, 50, 600, 100, 2, 10, 4)
    noisetest6 = (100, 80, 0, 40, 500, 20, 5, 10, 4)
    newTest = noisetest1
    
    noiseArgs = Arguments(newTest...)
    
    Smax = noiseArgs.maxStock
    T = noiseArgs.horizon
    
    dh_bellman_values = zeros(Smax+1, T+1)
    hd_bellman_values = zeros(Smax+1, T+1)
    hd_bellman_values_2 = zeros(Smax+1, T+1);
    
    xaxis, y1, y2 = olddependance_noise!(noiseArgs, dh_bellman_values, hd_bellman_values, hd_bellman_values_2)
    
    yaxis = [y1, y2]
    
    gr()
    
    display( plot(
    xaxis, yaxis, label=["case 1" "case 2"], legend = :bottomleft, legendfontsize=10,
    line=[:auto :auto],
    color=["blue" "red"], alpha=[0.9 0.9], xlabel="Max Noise "*L"\overline{W}"*" "*L"T="*"$horizon", xlabelfontsize=16,
    ylabel="Relative distance", ylabelfontsize=16,
    legendtitle="Cases", legendtitlefontsize=12,
    linewidth=3, thickness_scaling = 1, framestyle = :origin
    )
    )
    
    savefig("img/dependencenoise-$(newTest...).png")
end

function dependance_noise_test(horizon::Int64)
    noisetest1 = (100, 80, 0, 20, 500, 20, 2, 10, horizon)
    noisetest2 = (200, 160, 0, 20, 500, 20, 6, 25, 4)
    noisetest3 = (50, 45, 0, 20, 500, 20, 10, 5, 4)
    noisetest4 = (150, 130, 0, 20, 500, 20, 2, 10, 4)
    noisetest5 = (150, 130, 0, 50, 600, 100, 2, 10, 4)
    noisetest6 = (100, 80, 0, 40, 500, 20, 5, 10, 4)
    newTest = noisetest1
    
    noiseArgs = Arguments(newTest...)
    
    Smax = noiseArgs.maxStock
    T = noiseArgs.horizon
    
    dh_bellman_values = zeros(Smax+1, T+1)
    hd_bellman_values = zeros(Smax+1, T+1)
    
    xaxis, y1, y2 = dependance_noise!(noiseArgs, dh_bellman_values, hd_bellman_values)
    
    yaxis = [y1, y2]
    
    gr()
    
    display( plot(
    xaxis, yaxis, label=["DH" "HD"], legend = :bottomleft, legendfontsize=10,
    line=[:auto :auto],
    color=["blue" "red"], alpha=[0.9 0.9], xlabel="Max Noise "*L"\overline{W}"*" for "*L"T="*"$horizon", xlabelfontsize=16,
    ylabel=L"V_0(s)"*" for "*L"s = \overline{S}/2", ylabelfontsize=16,
    legendtitle="Cases", legendtitlefontsize=12,
    linewidth=3, thickness_scaling = 1, framestyle = :origin
    )
    )
    
    savefig("img/valuedepnoise$(newTest...).png")
end

#multiple_simulations_twodiff()
#oldmultiple_simulations()
#main()
@time olddependance_noise_test(parse(Int64, ARGS[1]))
#showPolicies()
