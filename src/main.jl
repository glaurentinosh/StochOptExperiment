include("dynProg.jl")
include("plotting.jl")
using Profile

function instantaneousCost(u, s, w, p)
    if (u <= s+w)
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

function dependancy_horizon(T, Smax, dh_bellman_values, hd_bellman_values)
    xaxis = 0:T
    #y1 = view(dh_bellman_values, round(Int, Smax/2), :)
    y1 = [i==T+1 ? 0 : dh_bellman_values[round(Int, Smax/2), i]/(T-i+1) for i in T+1:-1:1]
    #y2 = view(hd_bellman_values, round(Int, Smax/2), :)
    y2 = [i==T+1 ? 0 : hd_bellman_values[round(Int, Smax/2), i]/(T-i+1) for i in T+1:-1:1]
    yaxis = [y1, y2]
    return xaxis, yaxis
end

function dependancy_horizon2(T, Smax, dh_bellman_values, hd_bellman_values)
    xaxis = 0:T
    #y1 = view(dh_bellman_values, round(Int, Smax/2), :)
    y1 = [i==T+1 ? 0 : dh_bellman_values[round(Int, Smax/2), i]/(T-i+1) for i in T+1:-1:1]
    #y2 = view(hd_bellman_values, round(Int, Smax/2), :)
    y2 = [i==T+1 ? 0 : hd_bellman_values[round(Int, Smax/2), i]/(T-i+1) for i in T+1:-1:1]
    yaxis = [100*(y2[i] - y1[i])/y2[i] for i in 1:length(y1)] 
    return xaxis, yaxis
end

function mainhorizon(args, dhbellman, hdbellman)
    fillvalues!(args, dhbellman)
    fillvalues!(args, hdbellman)

    Smax = args.maxStock
    T = args.horizon

    xaxis, yaxis = dependancy_horizon(T, Smax, dhbellman, hdbellman)

    gr()
    
    display( plot(
    xaxis, yaxis, label=["DH" "HD"], legend = :bottomleft, legendfontsize=10,
    line=[:auto :auto],
    color=["blue" "red"], alpha=[0.9 0.9], xlabel="Horizon "*L"T", xlabelfontsize=16,
    ylabel=L"V_0(s)"*" for "*L"s = \overline{S}/2", ylabelfontsize=16,
    legendtitle="Cases", legendtitlefontsize=12,
    linewidth=3, thickness_scaling = 1, framestyle = :origin
    )
    )
    
    savefig("img/valuedephorizon$(name(args)).png")
end

function oldsimulate2(args::Arguments)
    Smax = args.maxStock
    
    dhbellman = BellmanFunctions{DH}(args)
    hdbellman = BellmanFunctions{HD}(args)
    
    fillvalues!(args, dhbellman)
    fillvalues!(args, hdbellman)
    
    return (hdbellman[round(Int, Smax/2),1] - dhbellman[round(Int, Smax/2),1])/hdbellman[round(Int, Smax/2),1]
end

function multipletests()
    diff = 0
    smax, umax, wmax, pmax, Pmax = 0,0,0,0,0
    umin = 0
    x = 0
    ustep = 20
    sstart = 50
    horizon = 8
    
    p = 10
    P = 0

    minNoise = 0
    multiplier = 0.25

    p = 10
    P = 0

    instCost(u, s, w) = instantaneousCost(u, s, w, p)
    finalCost(s) = finalCost_aux(s, P)
    
    for s in sstart:50:300, u in 0.5*s:ustep:0.9*s, w in 0.2*s:ustep:4*s#, p in 2:4:10#, P in 5:5:25
        println("s, u, w : $(s), $(u), $(w)")
        dynamics(x, u, w) = dynamics_aux(x, u, w, s)
        noise = Noise{NormalDistribution}(minNoise, round(Int, w), multiplier)
        args = Arguments(s,u,umin,1,noise, dynamics, instCost, finalCost, horizon)
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

function mainplots()
    maxStock = 50
    maxControl = 45
    minControl = 0
    stepControl = 1
    horizon = 8

    minNoise = 0
    maxNoise = 70
    multiplier = 0.25
    noise = Noise{NormalDistribution}(minNoise, maxNoise, multiplier)

    p = 10
    P = 0

    dynamics(s, u, w) = dynamics_aux(s, u, w, maxStock)
    instCost(u, s, w) = instantaneousCost(u, s, w, p)
    finalCost(s) = finalCost_aux(s, P)

    args = Arguments(maxStock, maxControl, minControl, stepControl, noise, dynamics, instCost, finalCost, horizon)

    dhbellman = BellmanFunctions{DH}(args)
    hdbellman = BellmanFunctions{HD}(args)

    #mainplot(args, dhbellman, hdbellman)
    maindepnoise(args, dhbellman, hdbellman, multiplier)
    #mainhorizon(args, dhbellman, hdbellman)
end

function main()
    mainplots()
    #multipletests()    
end

function benchmark()
    println("======================= First run:")
    @time main()
 
    println("\n\n======================= Second run:")
    Profile.init(delay=0.1)
    Profile.clear()
    Profile.clear_malloc_data()
    @profile @time main()
 
    r = Profile.retrieve()
    f = open("output_file_name.txt", "w")
    #Profile.print(f; combine=true, mincount=100)
    Profile.print(IOContext(f, :displaysize => (200,500)); mincount=100,combine=true)
    close(f)
 end

 benchmark()
