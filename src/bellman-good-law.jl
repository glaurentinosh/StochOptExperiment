using Plots
using LaTeXStrings
using Printf
using Distributions

const FloatInt = Union{Float64, Int64}

const instCostConst = 5
const finalCostConst = 10

const MAXSTOCK = 50
const MAXCONTROL = 40

mutable struct Variables
    maxStock::Int64
    maxControl::Int64
    minControl::Int64
    maxNoise::Int64
    horizon::Int64
    law::Array{FloatInt}
    instCost::Function
    finalCost::Function
    dynamics::Function
end

mutable struct NoiseArgs
    variables::Variables
    minNoise::Int64
    maxNoise::Int64
    stepNoise::Int64
end

function filename(vars::Variables)
    "$(vars.maxStock)-"*
    "$(vars.maxControl)-"*
    "$(vars.minControl)-"*
    "$(vars.maxNoise)-"*
    "$(vars.horizon)-"*
    "$(vars.instCost(1, 1, 1))-"*
    "$(vars.finalCost(1))"
end

function filename(noiseargs::NoiseArgs)
    varname = filename(noiseargs.variables)
    "$(noiseargs.minNoise), $(noiseargs.maxNoise), $(noiseargs.stepNoise);$(varname)"
end

function instCoststandard(control, stock, noise)
    if control <= stock
        return instCostConst*(control)
    else
        return -Inf
    end
end

function instCostnormalcenter(u, s, w)
    u <= s+w ? instCostConst*pdf(Normal(MAXCONTROL/2, MAXCONTROL/8), u) : -Inf
end

instCost(control, stock, noise) = instCostnormalcenter(control, stock, noise)

finalcostsine(s) = sin(pi*s/100)
finalCostlinear(s) = finalCostConst*s
finalcostnormal(s) = pdf(Normal(MAXSTOCK, MAXSTOCK/8), s)
finalcostnegative(s) = finalCostConst*(MAXSTOCK - s)
finalcostcossine(s) = finalCostConst*(1+0.5cos(2pi*s/MAXSTOCK))
finalcostnormalcenter(s) = pdf(Normal(MAXSTOCK/2, MAXSTOCK/8), s)
finalCost(s) = finalcostnormalcenter(s)

function fillvaluesdh!(vars::Variables, dhvaluefunction)
    maxStock = vars.maxStock
    maxControl = vars.maxControl
    minControl = vars.minControl
    maxNoise = vars.maxNoise
    law = vars.law
    horizon = vars.horizon
    instCost = vars.instCost
    finalCost = vars.finalCost
    dynamics = vars.dynamics
    
    for stock in 0:maxStock
        dhvaluefunction[stock+1,horizon+1] = finalCost(stock)
    end

    for t in horizon-1:-1:0, stock in 0:maxStock
        maxexpected = -Inf
        for control in minControl:maxControl
            expectedvalue = 0
            divideby = 0
            for w in 0:maxNoise
                realization = instCost(control, stock, w) + dhvaluefunction[dynamics(stock, control, w) + 1, t+2]
                expectedvalue += realization * law[w+1,t+1]
                divideby += law[w+1,t+1]
            end
            expectedvalue /= divideby
            if maxexpected < expectedvalue
                maxexpected = expectedvalue
            end
        end
        dhvaluefunction[stock+1,t+1] = maxexpected
    end
end

function fillvalueshd!(vars::Variables, hdvaluefunction)
    maxStock = vars.maxStock
    maxControl = vars.maxControl
    minControl = vars.minControl
    maxNoise = vars.maxNoise
    law = vars.law
    horizon = vars.horizon
    instCost = vars.instCost
    finalCost = vars.finalCost
    dynamics = vars.dynamics
    
    for stock in 0:maxStock
        hdvaluefunction[stock+1,horizon+1] = finalCost(stock)
    end
    
    for t in horizon-1:-1:0, stock in 0:maxStock
        expectedmax = 0
        divideby = 0
        for w in 0:maxNoise
            maxvalue = -Inf
            for control in minControl:maxControl
                realization = instCost(control, stock, w) + hdvaluefunction[dynamics(stock, control, w) + 1, t+2]
                if maxvalue < realization
                    maxvalue = realization
                end
            end
            expectedmax += maxvalue * law[w+1,t+1]
            divideby += law[w+1,t+1]
        end
        expectedmax /= divideby
        hdvaluefunction[stock+1,t+1] = expectedmax 
    end 
end

dynamics(s, u, w) = round(Int, max(min(s - u + w, maxStock), 0))

uniflaw(maxNoise, horizon) = [1 for i in 0:maxNoise, j in 0:horizon]
linearlaw(maxNoise, horizon) = [i for i in 0:maxNoise, j in 0:horizon]
normallaw(maxNoise, horizon) = [pdf(Normal(maxNoise/2, maxNoise/2^j), i) for i in 0:maxNoise, j in 0:horizon]
quadraticlaw(maxNoise, horizon) = [i^2 for i in 0:maxNoise, j in 0:horizon]
paraboliclaw(maxNoise, horizon) = [i*(maxNoise - i) for i in 0:maxNoise, j in 0:horizon]
parabolic2law(maxNoise, horizon) = [(i - maxNoise/2)^2 for i in 0:maxNoise, j in 0:horizon]
timedependentlaw(maxNoise, horizon) = [(1.5+cos(4pi*j/horizon))*(i-maxNoise/2)^2 for i in 0:maxNoise, j in 0:horizon]
randomlaw(maxNoise, horizon) = rand(maxNoise+1, horizon+1)
cossinelaw(maxNoise, horizon) = [1+cos(2pi*i/maxNoise) for i in 0:maxNoise, j in 0:horizon]
sinelaw(maxNoise, horizon) = [1+sin(2pi*i/maxNoise) for i in 0:maxNoise, j in 0:horizon]

normalizelaw(law) = law./sum(law, dims=1)

function plotlaw(thislaw)
    law = normalizelaw(thislaw)
    x_axis = 0:size(law)[1]-1
    y_axis = [view(law, :, i) for i in 1:size(law)[2]]
    
    display( plot(
    x_axis, y_axis, legend = :bottomright, legendfontsize=14,
    line=[:auto :auto], linewidht=[20 4],
    alpha=[0.9 0.9], xlabel="Noise "*L"(w)", xlabelfontsize=16,
    ylabel="Normalized cost " * latexstring("V_{t}(s)"), ylabelfontsize=16,
    legendtitle=latexstring("t = {t}"), legendtitlefontsize=16,
    linewidth=3, thickness_scaling = 1, framestyle = :origin
    #title="Normalized cost at time t = $t", titlefontsize=18
    )
    )
    
    #savefig("value_functions_T-$(T-t).png")
end

function normalize(vec, N)
    if N != 0
        return vec/N
    end
    return vec
end;

function plotting(T, Smax,dh_bellman,hd_bellman)
    for t in T:-1:0
        x_axis = 0:Smax
        y_axis = [
        normalize(view(dh_bellman, :, t+1), (T-t)),
        normalize(view(hd_bellman, :, t+1), (T-t)),
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

function noisedependancy!(
    args::NoiseArgs,dhbellman::Array{U, 2},hdbellman::Array{U, 2}
    ) where U <: FloatInt
    
    y1 = FloatInt[]
    y2 = FloatInt[]
    
    Smax = args.variables.maxStock
    x_axis = args.minNoise:args.stepNoise:args.maxNoise
    for w in x_axis
        args.variables.maxNoise = w
        println(w)
        
        @time fillvaluesdh!(args.variables, dhbellman)
        @time fillvalueshd!(args.variables, hdbellman)
        
        hd_middle_value = hdbellman[round(Int, Smax/2),1]
        dh_middle_value = dhbellman[round(Int, Smax/2),1]
    
        push!(y1, dh_middle_value)
        push!(y2, hd_middle_value)
    end
    
    return x_axis, y1, y2
end

function maindepnoise()
    thismaxStock = 50
    thismaxControl = 40
    thisminControl = 0
    thismaxNoise = 4*thismaxStock
    thishorizon = 8
    law = timedependentlaw(thismaxNoise, thishorizon)
    thislaw = normalizelaw(law)
    
    thisinstCost = instCost
    thisfinalCost = finalCost

    dynamics(s, u, w) = round(Int, max(min(s - u + w, thismaxStock), 0))
    thisdyn = dynamics

    vars = Variables(thismaxStock, thismaxControl, thisminControl, thismaxNoise,
                        thishorizon, thislaw, thisinstCost, thisfinalCost, thisdyn)

    dhvaluefunction = zeros(maxStock + 1, horizon + 1)
    hdvaluefunction = zeros(maxStock + 1, horizon + 1)
    
    noiseArgs = NoiseArgs(vars,0,thismaxNoise,20)
    
    xaxis, y1, y2 = noisedependancy!(noiseArgs, dhvaluefunction, hdvaluefunction)
    
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
    savefig("img/valuedepnoise$(filename(noiseArgs)).png")
end

function maindephorizon()
    thismaxStock = 50
    thismaxControl = 40
    thisminControl = 0
    thismaxNoise = 4*thismaxStock
    thishorizon = 50
    law = uniflaw(thismaxNoise, thishorizon)
    thislaw = normalizelaw(law)
    thisinstCost = instCost
    thisfinalCost = finalCost
    dynamics(s, u, w) = round(Int, max(min(s - u + w, thismaxStock), 0))
    thisdyn = dynamics

    vars = Variables(thismaxStock, thismaxControl, thisminControl, thismaxNoise,
                        thishorizon, thislaw, thisinstCost, thisfinalCost, thisdyn)

    dhvaluefunction = zeros(thismaxStock + 1, thishorizon + 1)
    hdvaluefunction = zeros(thismaxStock + 1, thishorizon + 1)

    fillvaluesdh!(vars, dhvaluefunction)
    fillvalueshd!(vars, hdvaluefunction)

    #xaxis, yaxis = dependancy_horizon(T, Smax, dh_bellman_values, hd_bellman_values)
    xaxis, yaxis = dependancy_horizon2(thishorizon, thismaxStock, dhvaluefunction, hdvaluefunction)

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
    
    savefig("img/valuedephorizon$(filename(vars)).png")
end

function mainplot(vars::Variables)  
    thismaxStock = vars.maxStock
    thishorizon = vars.horizon

    dhvaluefunction = zeros(thismaxStock + 1, thishorizon + 1)
    hdvaluefunction = zeros(thismaxStock + 1, thishorizon + 1)
    
    fillvaluesdh!(vars, dhvaluefunction)
    fillvalueshd!(vars, hdvaluefunction)

    plotting(thishorizon, thismaxStock, dhvaluefunction, hdvaluefunction)
end

function main()
    thismaxStock = 50
    thismaxControl = 40
    thisminControl = 0
    thismaxNoise = 40
    thishorizon = 50
    
    law = uniflaw(thismaxNoise,thishorizon)
    thislaw = normalizelaw(law)
    
    thisinstCost = instCost
    thisfinalCost = finalCost

    dynamics(s, u, w) = round(Int, max(min(s - u + w, thismaxStock), 0))
    thisdyn = dynamics
    
    vars = Variables(thismaxStock, thismaxControl, thisminControl, thismaxNoise,
    thishorizon, thislaw, thisinstCost, thisfinalCost, thisdyn)

    #maindepnoise()
    mainplot(vars)
    #maindephorizon()
    #plotlaw(sinelaw(40, 4))

end

main()

