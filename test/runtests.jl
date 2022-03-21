using TunabilityOscillations
using DifferentialEquations
using Test

using Plots

#=
I think I should define a function for each test to make it more modular and
easier to understand. Right now the code is a little messy because it has all
the code within the testset. I should also name the test function so it
reveals what is being tested. I'm not sure if the test function should return
a value or the logical test. Maybe returning true of false would be OK. Since
someone reading the code would go to that function to see what is really
happening
=#

@testset "Model creation" begin include("test_model_creation.jl") end

#=
## Test that the created model equations behave like models defined by hand

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
