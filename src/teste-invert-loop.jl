using BenchmarkTools

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

function Lt(u, s, w, args::Arguments)
    Smax = args.maxStock
    p = args.instantaneousCostConstant
    if (u <= min(s+w, Smax))
        return p*u
    else
        return -Inf
    end
end

function dynamics_aux(s, u, w, args::Arguments)
    Smax = args.maxStock
    return round(Int, max(min(s - u + w, Smax), 0))
end

function fill_bad(hd_bellman_function::Array{Float64, 2}, args::Arguments)
    L(u, s, w) = Lt(u, s, w, args)
    dynamics(s, u, w) = dynamics_aux(s, u, w, args)

    Smax = args.maxStock
    Umax = args.maxControl
    Wmax = args.maxNoise
    p = args.instantaneousCostConstant
    T = args.horizon
    Wmin = 0
    
    # Calculate hd bellman functions    
    for t in T-1: -1: 0 # t = T-1, ..., 0
        w = 0 # avoid creating local variables
        for s in 0:Smax # s = 0, 1, ..., Smax
            expected_value = 0
            for i in 0:Wmax # i = 0, 1, ..., Nw
                w = Wmin + i * (Wmax - Wmin)
                max_value = -Inf
                for u in 0:Umax # u = Umin, ..., Umax 
                    max_value = max(max_value, L(u, s, w) + hd_bellman_function[t+2,dynamics(s, u, w)+1])
                end
                expected_value += max_value
            end
            expected_value /= (Wmax + 1)
            hd_bellman_function[t+1,s+1] = expected_value
        end
    end
end

function fill_good(hd_bellman_function::Array{Float64, 2}, args::Arguments)
    L(u, s, w) = Lt(u, s, w, args)
    dynamics(s, u, w) = dynamics_aux(s, u, w, args)

    Smax = args.maxStock
    Umax = args.maxControl
    Wmax = args.maxNoise
    p = args.instantaneousCostConstant
    T = args.horizon
    Wmin = 0
    
    # Calculate hd bellman functions    
    for t in T-1: -1: 0 # t = T-1, ..., 0
        for s in 0:Smax # s = 0, 1, ..., Smax
            expected_value = 0
            for i in 0:Wmax # i = 0, 1, ..., Nw
                w = Wmin + i * (Wmax - Wmin)
                max_value = -Inf
                for u in 0:Umax # u = Umin, ..., Umax 
                    max_value = max(max_value, L(u, s, w) + hd_bellman_function[dynamics(s, u, w)+1, t+1])
                end
                expected_value += max_value
            end
            expected_value /= (Wmax + 1)
            hd_bellman_function[s+1,t+1] = expected_value
        end
    end
end

function fill_good_dontavoidlocals(hd_bellman_function::Array{Float64, 2}, args::Arguments)
    L(u, s, w) = Lt(u, s, w, args)
    dynamics(s, u, w) = dynamics_aux(s, u, w, args)

    Smax = args.maxStock
    Umax = args.maxControl
    Wmax = args.maxNoise
    p = args.instantaneousCostConstant
    T = args.horizon
    Wmin = 0
    
    # Calculate hd bellman functions    
    for t in T-1: -1: 0 # t = T-1, ..., 0
        for s in 0:Smax # s = 0, 1, ..., Smax
            expected_value = 0
            for i in 0:Wmax # i = 0, 1, ..., Nw
                w = Wmin + i * (Wmax - Wmin)
                max_value = -Inf
                for u in 0:Umax # u = Umin, ..., Umax 
                    @inbounds begin max_value = max(max_value, L(u, s, w) + hd_bellman_function[dynamics(s, u, w)+1, t+1]) end
                end
                expected_value += max_value
            end
            expected_value /= (Wmax + 1)
            hd_bellman_function[s+1,t+1] = expected_value
        end
    end
end

function main()
    args = Arguments(2000,2,0,1,2,0,5,10,2000)
    hd_bad = zeros(args.horizon + 1, args.maxStock + 1)
    hd_good = zeros(args.maxStock + 1, args.horizon + 1)
    hd_good_dontavoidlocals = zeros(args.maxStock + 1, args.horizon + 1)

    t = @elapsed fill_bad(hd_bad, args)
    println(rpad("laplacian_bad", 30), t)
    
    t = @elapsed fill_good(hd_good, args)
    println(rpad("laplacian_good", 30), t)

    t = @elapsed fill_good_dontavoidlocals(hd_good_dontavoidlocals, args)
    println(rpad("laplacian_good_dntavdlcls", 30), t)
end
