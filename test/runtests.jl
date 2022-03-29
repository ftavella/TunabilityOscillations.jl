using TunabilityOscillations
# using DifferentialEquations
using Catalyst, LSODA
using Test
using Plots

@testset "Model creation" begin include("test_model_creation.jl") end

#=
# Goodwin Oscillator

# Test that the period of the oscillations is correct
# TODO
tspan = (0.0, 1.0)
prob = ODEProblem(model_equations!, u0, tspan, p)
ode_solution = solve(prob, isoutofdomain=(y,p,t)->any(x->x<0,y))

f_sampling = 4000
t_eval = LinRange(0.0, 1.0, f_sampling)
pks = find_periodogram_peak(ode_solution, t_eval, f_sampling)
# println(pks)
=#
