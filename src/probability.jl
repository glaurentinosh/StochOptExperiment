using Distributions

abstract type ProbabilityDistribution end

abstract type LinearDistribution <: ProbabilityDistribution end
abstract type UniformDistribution <: ProbabilityDistribution end
abstract type QuadraticDistribution <: ProbabilityDistribution end
abstract type Quadratic2Distribution <: ProbabilityDistribution end
abstract type NormalDistribution <: ProbabilityDistribution end

abstract type RandomVariable{T<:ProbabilityDistribution} end

mutable struct Noise{T} <: RandomVariable{T}
    noiseRange::UnitRange{Int64}
    law::Function
end

function Noise(oldNoise::Noise{T}, newMaxNoise::Int64, newmultiplier = 1) where T<:ProbabilityDistribution
    Noise{T}(first(oldNoise.noiseRange), newMaxNoise, newmultiplier)
end

function Noise{LinearDistribution}(minNoise::Int64, maxNoise::Int64, multiplier = 1)
    noiseRange = minNoise:maxNoise
    law(w::Int64) = linearLaw(w, minNoise, maxNoise, multiplier)
    Noise{LinearDistribution}(noiseRange, law)
end

function Noise{UniformDistribution}(minNoise::Int64, maxNoise::Int64, multiplier = 1)
    noiseRange = minNoise:maxNoise
    law(w::Int64) = uniformLaw(minNoise, maxNoise)
    Noise{UniformDistribution}(noiseRange, law)
end

function Noise{QuadraticDistribution}(minNoise::Int64, maxNoise::Int64, multiplier = 1)
    noiseRange = minNoise:maxNoise
    law(w::Int64) = quadraticLaw(w, minNoise, maxNoise)
    Noise{QuadraticDistribution}(noiseRange, law)
end

function Noise{Quadratic2Distribution}(minNoise::Int64, maxNoise::Int64, multiplier = 1)
    noiseRange = minNoise:maxNoise
    law(w::Int64) = quadratic2Law(w, minNoise, maxNoise)
    Noise{Quadratic2Distribution}(noiseRange, law)
end

function Noise{NormalDistribution}(minNoise::Int64, maxNoise::Int64, multiplier = 1)
    noiseRange = minNoise:maxNoise
    law(w::Int64) = normalLaw(w, minNoise, maxNoise, multiplier)
    Noise{NormalDistribution}(noiseRange, law)
end

function linearLaw(w::Int64, minNoise::Int64, maxNoise::Int64, multiplier::Real)
    minProb = 2/((multiplier+1)*(maxNoise - minNoise + 1))
    deltaProb = (multiplier - 1) * minProb / (maxNoise - minNoise)
    return minProb + deltaProb * (w - minNoise)
end

function uniformLaw(minNoise::Int64, maxNoise::Int64)
    return 1/(maxNoise - minNoise + 1)
end

function quadraticLaw(w::Int64, minNoise::Int64, maxNoise::Int64)
    return 6*(w - minNoise)*(maxNoise - w)/(maxNoise - minNoise)^3
end

function quadratic2Law(w::Int64, minNoise::Int64, maxNoise::Int64)
    return 12(w - (minNoise + maxNoise)/2)^2/(maxNoise - minNoise)^3
end


function normalLaw(w::Int64, minNoise::Int64, maxNoise::Int64, multiplier::Real)
    mean = (maxNoise - minNoise)/2
    std = multiplier * (maxNoise - mean)
    n = Normal(mean, std)
    return pdf(n, w)
end