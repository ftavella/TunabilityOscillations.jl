"""
    calc_main_freq(ode_sol::ODESolution, sol_sampling::Int)

Compute the main Fourier frequency and its power for each node of an ODE solution.

Input sol_sampling determines how many equally-spaced points are used to interpolate
the ODE solution.
"""
function calc_main_freq(ode_sol::ODESolution, sol_sampling::Int)
  N = length(ode_sol.u[1])
  freq = zeros(N)
  power = zeros(N)
  t_eval = LinRange(ode_sol.t[1], ode_sol.t[end], sol_sampling)
  fs = round(Int, sol_sampling/ode_sol.t[end])
  for i in 1:N
    # Remove mean to avoid peak at low frequencies
    signal = ode_sol(t_eval)[i,:] .- mean(ode_sol(t_eval)[i,:])
    pgram = periodogram(signal; fs=fs)
    pks, ~ = findmaxima(pgram.power)
    pks, ~ = peakproms(pks, pgram.power, minprom=maximum(pgram.power)/3.0)
    if length(pks) > 1
      # From the top peaks get the one with the lowest frequency
      best_vals = sort(collect(zip(pgram.freq[pks], pgram.power[pks])); by=first)
      freq[i] = best_vals[1][1]
      power[i] = best_vals[1][2]
    else
      freq[i] = pgram.freq[argmax(pgram.power)]
      power[i] = pgram.power[argmax(pgram.power)]
    end
  end
  return [freq, power]
end


"""
    calc_amp(ode_sol::ODESolution, peak_num_thresh::Int)

Compute the peak-to-peak distance for all variables in an ODE solution and the
coefficient of variation for maxima and minima peaks.

Input peak_num_thresh is an integer used to filter non-oscillatory solutions.
"""
function calc_amp(ode_sol::ODESolution, peak_num_thresh::Int)
  N = length(ode_sol.u[1])
  amp_val = zeros(N)
  amp_cv_max = zeros(N)
  amp_cv_min = zeros(N)
  for i in 1:N
    rough_amp = maximum(ode_sol[i,:]) - minimum(ode_sol[i,:])
    pks_max, ~ = findmaxima(ode_sol[i,:])
    pks_max, ~ = peakproms(pks_max, ode_sol[i,:], minprom=rough_amp/3.0)
    vals_max = ode_sol[i, pks_max]
    pks_min, ~ = findminima(ode_sol[i,:])
    pks_min, ~ = peakproms(pks_min, ode_sol[i,:], minprom=rough_amp/3.0)
    vals_min = ode_sol[i, pks_min]
    if length(vals_max) > peak_num_thresh
      if length(vals_min) > peak_num_thresh
        if length(vals_max) < 3*peak_num_thresh
          if length(vals_min) < 3*peak_num_thresh
            amp_val[i] = maximum(vals_max) - minimum(vals_min)
            amp_cv_max[i] = std(vals_max)/mean(vals_max)
            amp_cv_min[i] = std(vals_min)/mean(vals_min)
          end
        end
      end
    end
  end
  return [amp_val, amp_cv_max, amp_cv_min]
end

"""
    ode_params_from_param_set(model::ReactionSystem, p_set::Vector)

Reorder and normalize an array of parameters to be suitable for remaking an
ODEProblem.

The returned parameter array has the right order to use on ODE solvers.
Input variable p_set is assumed to be ordered as:
[α1, α2, ..., β1, ..., γ1, ..., κ1, ..., η1, ..., ηN]
"""
function ode_params_from_param_set(model::ReactionSystem, p_set::Vector)
  N = length(species(model))
  P = length(parameters(model))
  E = Int((P - 2N)/3)
  new_map_vals = []
  for j in 1:N
    @nonamespace push!(new_map_vals, (model.α[j], p_set[j]))
    @nonamespace push!(new_map_vals, (model.β[j], p_set[N+j]))
  end
  for j in 1:E
    @nonamespace push!(new_map_vals, (model.γ[j], p_set[2*N+j]))
    @nonamespace push!(new_map_vals, (model.κ[j], p_set[2*N+E+j]))
    @nonamespace push!(new_map_vals, (model.η[j], p_set[2*N+2*E+j]))
  end
  return ModelingToolkit.varmap_to_vars(Dict(new_map_vals), parameters(model))
end


"""
    generate_LHC_sample(model::ReactionSystem, samples::Int, lim::Dict)

Create an array of parameter sets by LatinHypercube or Random sampling.

A total of `samples` parameter sets are generated. Sample values are scaled by
the quantities defined in `lim`. Scaling can be linear or logarithmic for α, β, and γ.
"""
function generate_LHC_sample(model::ReactionSystem, samples::Int, lim::Dict, hparams::Dict)
  N = length(species(model))
  P = length(parameters(model))
  E = Int((P - 2N)/3)
  # Setup scaling rule for the sample
  scaling_plan = Tuple{Float64,Float64}[]
  if lowercase(hparams["sampling_scale"]) == "log"
    scaling_plan = vcat(scaling_plan, [(log10(lim["α"][1]), log10(lim["α"][2])) for i in 1:N])
    scaling_plan = vcat(scaling_plan, [(log10(lim["β"][1]), log10(lim["β"][2])) for i in 1:N])
    scaling_plan = vcat(scaling_plan, [(log10(lim["γ"][1]), log10(lim["γ"][2])) for i in 1:E])
  elseif lowercase(hparams["sampling_scale"]) == "linear"
    scaling_plan = vcat(scaling_plan, [(lim["α"][1], lim["α"][2]) for i in 1:N])
    scaling_plan = vcat(scaling_plan, [(lim["β"][1], lim["β"][2]) for i in 1:N])
    scaling_plan = vcat(scaling_plan, [(lim["γ"][1], lim["γ"][2]) for i in 1:E])
  end
  scaling_plan = vcat(scaling_plan, [(lim["κ"][1], lim["κ"][2]) for i in 1:E])
  scaling_plan = vcat(scaling_plan, [(lim["η"][1], lim["η"][2]) for i in 1:E])

  # Draw sample
  if lowercase(hparams["sampling_style"]) == "lhc"
    LHCplan = randomLHC(samples, P)
  elseif lowercase(hparams["sampling_style"]) == "random"
    LHCplan = rand([i for i in 1:samples], samples, P)
  end

  LHC = scaleLHC(LHCplan, scaling_plan)
  # Rescale accordingly if sampling scale is log
  if lowercase(hparams["sampling_scale"]) == "log"
    LHC[:, 1:2*N+E] = 10 .^ LHC[:, 1:2*N+E]
  end

  # Make it dimensionless in time
  LHC[:, 1] .= 1.0

  return LHC
end

"""
    create_generic_ode_problem(model::ReactionSystem)

Create a generic ODE problem from a model.
"""
function create_generic_ode_problem(model::ReactionSystem)
  N = length(species(model))
  P = length(parameters(model))
  E = Int((P - 2N)/3)
  p_dict_vals = []
  @nonamespace p_dict_vals = vcat(p_dict_vals, [(model.α[i], 0.5) for i in 1:N])
  @nonamespace p_dict_vals = vcat(p_dict_vals, [(model.β[i], 1.0) for i in 1:N])
  @nonamespace p_dict_vals = vcat(p_dict_vals, [(model.γ[i], 1.0) for i in 1:E])
  @nonamespace p_dict_vals = vcat(p_dict_vals, [(model.κ[i], 0.5) for i in 1:E])
  @nonamespace p_dict_vals = vcat(p_dict_vals, [(model.η[i], 1.0) for i in 1:E])
  pmap = Dict(p_dict_vals)
  @nonamespace u0_dict_vals = [(model.a[i], 0.5) for i in 1:N]
  u₀map = Dict(u0_dict_vals)
  return ODEProblem(model, u₀map, (0.0, 1.0), pmap)
end

"""
    infer_timescale(model::ReactionSystem, p_set::Vector)

Calculate the longest timescale of a model as the minimum of the reaction rates.
"""
function infer_timescale(model::ReactionSystem, p_set::Vector)
  N = length(species(model))
  P = length(parameters(model))
  E = Int((P - 2N)/3)
  # Only α through γ affect the timescale of the solution.
  rates = [p_set[k] for k in 1:(2*N+E)]
  timescale = 1.0/minimum(rates)
  return timescale
end

"""
    is_oscil(freq::Vector, power::Vector, amp::Vector, amp_cv_max::Vector,
             amp_cv_min::Vector, hparams::Dict)

Decide if the features of an ODE solution indicate oscillations.
"""
function is_oscil(freq::Vector, power::Vector, amp::Vector, amp_cv_max::Vector,
                  amp_cv_min::Vector, hparams::Dict)
    N = length(freq)
    oscillatory = false
    freqtol = hparams["freq_tolerance"]
    freq_diffs = [abs(freq[j] - freq[j+1])/freq[j] < freqtol for j in 1:N-1]
    if all(freq_diffs)
      if all([power[j] > hparams["power_threshold"] for j in 1:N])
        if any(amp_cv_max .< hparams["amp_cv_thresh"])
          if any(amp_cv_min .< hparams["amp_cv_thresh"])
            if all(amp .> 0.0)
              oscillatory = true
            end
          end
        end
      end
    end
    return oscillatory
end

"""
    find_oscillations(model::ReactionSystem, samples::Int, param_limits::Dict,
                      hparams::Dict)

Compute parameter sets that produce oscillations in a model.
"""
function find_oscillations(model::ReactionSystem, samples::Int,
                           param_limits::Dict, hparams::Dict)
  N = length(species(model))
  p_sample = generate_LHC_sample(model, samples, param_limits, hparams)
  prob = create_generic_ode_problem(model)

  """
      prob_func(prob::ODEProblem, i, ~)

  Equilibrate parameter set and define both simulation time and starting point
  for the model's simulation using EnsembleProblem.
  """
  function prob_func(prob::ODEProblem, i, ~)
    sample_p_set = p_sample[Int(i),:]
    ode_p_set = ode_params_from_param_set(model, sample_p_set)
    #  Equilibrate parameter set
    timescale = infer_timescale(model, sample_p_set)
    t_final_eq = hparams["equil_tscales"]*timescale
    equil_prob = remake(prob, p=ode_p_set, u0=ones(N)/2.0,
                        tspan=[0.0, t_final_eq])
    equil_sol = solve(equil_prob, hparams["solver"], abstol=hparams["abstol"],
                      reltol=hparams["reltol"], maxiters=hparams["maxiters"],
                      dense=false)
    # Set parameters for simulation
    midpoint = round(Int, length(equil_sol.t)/2.0)
    final_half_eq = equil_sol[midpoint:end]
    sol_sampling = length(final_half_eq.t)
    freq, ~ = calc_main_freq(final_half_eq, sol_sampling)
    if minimum(freq) > 0.0
      new_timescale = 0.5/minimum(freq) # Because I only considered half the sim
    else
      new_timescale = timescale
    end
    t_final = hparams["sim_tscales"] * new_timescale
    remake(prob, p=ode_p_set, u0=equil_sol.u[end], tspan=[0.0, t_final])
  end


  """
      output_func(sol::ODESolution, ~)

  Calculate frequency and amplitude and determine if solution is oscillatory.
  """
  function output_func(sol::ODESolution, ~)
    sol_sampling = length(sol.t)
    freq, power = calc_main_freq(sol, sol_sampling)
    amp, amp_cv_max, amp_cv_min = calc_amp(sol, hparams["peak_num_thresh"])
    oscillatory = is_oscil(freq, power, amp, amp_cv_max, amp_cv_min, hparams)
    endpoint = sol.u[end]
    out = [oscillatory, endpoint, freq, amp]
    return (out, false)
  end

  # Setup EnsembleProblem and solve it
  ensemble_prob = EnsembleProblem(prob, prob_func=prob_func, safetycopy=false,
                                  output_func=output_func)
  sim = solve(ensemble_prob, hparams["solver"], EnsembleThreads(), 
              trajectories=samples, abstol=hparams["abstol"], reltol=hparams["reltol"], 
              maxiters=hparams["maxiters"], dense=false)
  # Return only data from oscillatory solutions
  osci_idxs = findall(x->x[1], sim.u)
  return [sim.u[osci_idxs], p_sample[osci_idxs,:]]
end
