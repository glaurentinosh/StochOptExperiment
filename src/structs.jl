include("probability.jl")

mutable struct Arguments
    maxStock::Int64
    maxControl::Int64
    minControl::Int64
    stepControl::Int64
    noise::Noise
    dynamics::Function
    instantaneousCost::Function
    finalCost::Function
    horizon::Int64
end

function Arguments(args::Arguments, newmaxnoise::Int64, multiplier = 1)
    newnoise = Noise(args.noise, newmaxnoise, multiplier)
    Arguments(args.maxStock, args.maxControl, 0, 1, newnoise, args.dynamics,
    args.instantaneousCost, args.finalCost, args.horizon)
end

function name(args::Arguments)
    return "$(args.maxStock)-$(args.maxControl)-$(last(args.noise.noiseRange))-"*
    "$(args.instantaneousCost(1,1,1))-$(args.finalCost(1))-$(args.horizon)"
end

abstract type ProblemType end

abstract type DH <: ProblemType end
abstract type HD <: ProblemType end

const FloatInt = Union{Float64, Int64}

abstract type ValueFunctions <: AbstractArray{FloatInt, 2} end

struct BellmanFunctions{T<:ProblemType} <: ValueFunctions
    data::Array{FloatInt, 2}
end

Base.getindex(bellman::BellmanFunctions, i::Int) = bellman.data[i]
Base.setindex!(bellman::BellmanFunctions, value, key::Int64) = (bellman.data[key] = value)
Base.size(bellman::BellmanFunctions) = size(bellman.data)
Base.IndexStyle(::Type{<:BellmanFunctions}) = IndexLinear()
Base.strides(bellman::BellmanFunctions) = strides(bellman.data)

function BellmanFunctions{T}(maxStock::Int64, horizon::Int64) where T <: ProblemType
    BellmanFunctions{T}(zeros(maxStock+1, horizon+1))
end

function BellmanFunctions{T}(args::Arguments) where T <: ProblemType
    BellmanFunctions{T}(zeros(args.maxStock + 1, args.horizon + 1))
end
