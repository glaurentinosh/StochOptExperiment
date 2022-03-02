struct Noise
    maxNoise::Int64
    minNoise::Int64
    law::Function
end

struct Arguments
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

abstract type ProblemType end

abstract type DH <: ProblemType end
abstract type HD <: ProblemType end

abstract type ValueFunctions <: AbstractArray{Real, 2} end

struct BellmanFunctions{T<:ProblemType} <: ValueFunctions
    data::Array{Real, 2}
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
