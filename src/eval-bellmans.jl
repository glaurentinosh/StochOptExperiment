# initialize bellman functions

function call(args...)
    Smax = args[1]
    Umax = args[2]
    Wmax = args[3]
    p = args[4]
    Pf = args[5]
    T = args[6]

    bellman_function, hd_bellman_function = args[7], args[8]
    bellman_function = []
    hd_bellman_function = []

#= 
    bellman_function = Array{Real}[]
    hd_bellman_function = Array{Real}[] =#

    # Possible states
    #Smax = parse(Int64, ARGS[1]) #200
    #Smin = 0
    
    # Possible decisions
    #Umax = parse(Int64, ARGS[2]) #80
    Umin = -Umax
    Ustep = 10
    
    # Possible scenarios
    #Wmax = parse(Int64, ARGS[3]) #50
    Wmin = 0
    ;
    
    # Instantaneous cost function
    #p = parse(Int64, ARGS[4]) #5
    
    function L(u, s, w)
        if (u <= min(s+w, Smax))
            return p*u#p*sqrt(u^2 + (s+w)^2)
        else
            return -Inf
        end
    end

    function L2(u, s, w)
        if (u <= s)
            return p*u#p*sqrt(u^2 + (s+w)^2)
        else
            return -Inf
        end
    end
    
    # Final cost
    #Pf = parse(Int64, ARGS[5]) #10
    function K1(s)
        if (s >= 0.8*Smax)
            return Pf*s
        else
            return 0
        end
    end
    
    function K(s)
        return Pf*s
    end
    
    # dynamics
    function dynamics(s, u, w)
        return round(Int, max(min(s - u + w, Smax), 0))
    end
    ;
    
    # horizon
    #T = parse(Int64, ARGS[6]) #8
    
    # number of scenarios
    Nw = Wmax
    Nw_inv = 1/Nw;
    
    for t in 1:T
        push!(bellman_function, zeros(Smax + 1))
        push!(hd_bellman_function, zeros(Smax + 1))
    end
    
    # final cost
    bellman_T = zeros(Smax + 1)
    hd_bellman_T = zeros(Smax + 1)
    
    for s in 1:Smax+1
        bellman_T[s] = K(s - 1)
        hd_bellman_T[s] = K(s - 1)
    end
    
    push!(bellman_function, bellman_T)
    push!(hd_bellman_function, hd_bellman_T)
    ;
    
    #println(typeof(bellman_function), size(bellman_function))
    #println(typeof(bellman_function[3]), size(bellman_function[3]))
    
    #bellman_function[T+1];
    
    # Calculate bellman functions    
    for t in T-1:-1:0 # t = T-1, ..., 0
        w = 0 #avoid creating local variables
        for s in 0:Smax # s = 0, 1, ..., Smax
            max_expected = -Inf # needs to be inf
            for u in Umin:Ustep:Umax # u = Umin, Umin+1, ..., Umax 
                expected_value = 0
                for i in 0:Nw # i = 0, 1, ..., Nw
                    w = Wmin + i * Nw_inv * (Wmax - Wmin)
                    expected_value += L(u, s, w) + bellman_function[t+2][dynamics(s, u, w) + 1]
                end
                expected_value /= (Nw + 1) # Nw + 1 realizations
                max_expected = max(expected_value, max_expected)
            end
            bellman_function[t+1][s+1] = max_expected
        end
    end
    
    
    # Calculate hd bellman functions    
    for t in T-1: -1: 0 # t = T-1, ..., 0
        w = 0 # avoid creating local variables
        for s in 0:Smax # s = 0, 1, ..., Smax
            expected_value = 0
            for i in 0:Nw # i = 0, 1, ..., Nw
                w = Wmin + i * Nw_inv * (Wmax - Wmin)
                max_value = -Inf
                for u in Umin:Ustep:Umax # u = Umin, ..., Umax 
                    max_value = max(max_value, L(u, s, w) + hd_bellman_function[t+2][dynamics(s, u, w)+1])
                end
                expected_value += max_value
            end
            expected_value /= (Nw + 1)
            hd_bellman_function[t+1][s+1] = expected_value
        end
    end

    return (hd_bellman_function[1][round(Int, Smax/2)] - bellman_function[1][round(Int, Smax/2)])/hd_bellman_function[1][round(Int, Smax/2)]

end