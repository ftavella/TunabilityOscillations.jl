module TunabilityOscillations

using Catalyst
using DifferentialEquations
using Statistics
using DSP
using Peaks
using LatinHypercubeSampling
using LSODA

export create_model, calc_main_freq, calc_amp, find_oscillations
export create_generic_ode_problem, ode_params_from_param_set

include("create_model_equations.jl")
include("detect_oscillations.jl")

end