module TunabilityOscillations

using Catalyst
using DifferentialEquations
using Statistics
using DSP
using Peaks
using LatinHypercubeSampling
using LSODA

export create_model, calc_main_freq, find_oscillations

include("create_model_equations.jl")
include("detect_oscillations.jl")

end
