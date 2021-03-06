using TunabilityOscillations
using Catalyst, LSODA
using Test
using Plots

@testset "Model creation" begin include("test_model_creation.jl") end
@testset "Frequency peak finding" begin include("test_periodogram_peak.jl") end
@testset "Oscillation detection" begin include("test_find_oscillations.jl") end
