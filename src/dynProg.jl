include("structs.jl")

function fillvalues!(args::Arguments, bellmanfunction::BellmanFunctions)
    maxStock = args.maxStock
    horizon = args.horizon
    finalCost(stock) = args.finalCost(stock)
    
    policies = zeros(horizon)

    for stock in 1:maxStock + 1
        bellmanfunction[stock, horizon + 1] = finalCost(stock - 1)
    end
    
    for t in horizon-1:-1:0, stock in 0:maxStock
        bellmanfunction[stock+1,t+1], policies[t+1] = dynamicProgramming!(t, stock, args, bellmanfunction)
    end

    return policies
end

function dynamicProgramming!(t::Int64, stock::Int64, args::Arguments, bellmanfunction::BellmanFunctions{DH})
    maxControl = args.maxControl
    minControl = args.minControl
    stepControl = args.stepControl
    noise = args.noise
    instCost(control, stock) = args.instantaneousCost(control, stock)
    dynamics(stock, control, w) = args.dynamics(stock, control, w)
    
    maxexpected = -Inf
    policy = 0

    for control in minControl:stepControl:maxControl
        expectedvalue = 0
        for w in noise.minNoise:noise.maxNoise
            realization = instCost(control, stock) + bellmanfunction[dynamics(stock, control, w) + 1, t+2]
            expectedvalue += realization * noise.law(w) 
        end
        if maxexpected < expectedvalue
            maxexpected = expectedvalue
            policy = control
        end
    end

    return maxexpected, policy
end

function dynamicProgramming!(t::Int64, stock::Int64, args::Arguments, bellmanfunction::BellmanFunctions{HD})
    maxControl = args.maxControl
    minControl = args.minControl
    stepControl = args.stepControl
    noise = args.noise
    instCost(control, stock) = args.instantaneousCost(control, stock)
    dynamics(stock, control, w) = args.dynamics(stock, control, w)
    
    expectedmax = 0
    policy = 0

    for w in noise.minNoise:noise.maxNoise
        maxvalue = -Inf
        for control in minControl:stepControl:maxControl
            realization = instCost(control, stock) + bellmanfunction[dynamics(stock, control, w) + 1, t+2]
            if maxvalue < realization
                maxvalue = realization
                policy = control
            end
        end
        expectedmax += maxvalue * noise.law(w)
    end

    return expectedmax, policy
end