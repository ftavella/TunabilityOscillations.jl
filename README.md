# Studying the Tunability of Oscillatory Networks

[![Build Status](https://github.com/ftavella/TunabilityOscillations.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/ftavella/TunabilityOscillations.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/ftavella/TunabilityOscillations.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ftavella/TunabilityOscillations.jl)

This package allows to:
- Define networks of interacting biomolecules through a connectivity matrix
- Sample parameter space using Latin Hypercube Sampling
- Discriminate between oscillatory and non-oscillatory solutions
- Calculate the frequency and amplitude of each oscillatory parameter set

For example, the following script calculates oscillatory solutions for the repressilator network
```julia
using TunabilityOscillations
using Catalyst, LSODA
using Plots, JLD2

samples = Int(1e6)
connectivity = [0 0 -1; -1 0 0; 0 -1 0]
hparams = Dict([("peak_num_thresh", 30), ("freq_tolerance", 0.01),
                ("equil_tscales", 50), ("sim_tscales", 40),
                ("power_threshold", 1e-7), ("amp_cv_thresh", 0.05),
                ("abstol", 1e-12), ("reltol", 1e-6)])
param_limits = Dict([
  ("α", [0.0, 1.0]),
  ("β", [10.0, 100.0]),
  ("γ", [100.0, 10000.0]),
  ("κ", [0.0, 1.0]),
  ("η", [3.0, 10.0]),])

model = create_model(connectivity)
@time sim, p_sample = find_oscillations(model, samples, param_limits, hparams)
jldsave("Repressilator_result.jld2"; samples, connectivity, hparams, param_limits, model, sim, p_sample)
```

Developed at the [Yang Lab - University of Michigan](http://www-personal.umich.edu/~qiongy/)
