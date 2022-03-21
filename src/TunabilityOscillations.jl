module TunabilityOscillations

# Write your package code here.
using DifferentialEquations
using Statistics
using DSP
using Peaks

export create_model, find_periodogram_peak

include("create_model_equations.jl")
include("detect_oscillations.jl")

end
