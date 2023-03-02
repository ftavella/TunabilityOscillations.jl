using TunabilityOscillations, Peaks
using DifferentialEquations, Statistics
using Test, Catalyst, Random, JLD2, Plots

Random.seed!(123)

# @testset "Model creation" begin include("test_model_creation.jl") end
# @testset "Frequency peak finding" begin include("test_periodogram_peak.jl") end
# @testset "Oscillation detection" begin include("test_find_oscillations.jl") end

@testset "Creating model from connectivity" begin include("test_create_model.jl") end
