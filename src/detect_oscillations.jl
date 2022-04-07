function find_periodogram_peak(ode_solution, t_eval, f_sampling)
  nodes = size(ode_solution, 1)
  peaks = fill(Float64[], nodes)
  for i in 1:nodes
    # Remove mean from solution to avoid peak at low frequencies
    signal = ode_solution(t_eval)[i,:] .- mean(ode_solution(t_eval)[i,:])
    pgram = periodogram(signal; fs=f_sampling)
    # Find and save largest peak in the periodogram
    pks, vals = findmaxima(pgram.power, 3)
    if isempty(pks) || isempty(vals)
      # If no peak is found, save maximum value of power spectrum
      peaks[i] = [pgram.freq[argmax(pgram.power)],
                  pgram.power[argmax(pgram.power)]]
    else
      peaks[i] = [pks[argmax(vals)], vals[argmax(vals)]]
    end
  end
  return peaks
end

#=
This function divides by β1 to create a dimensionless system
=#
function select_params_from_sample(model, p_sample, i)
  N = length(species(model))
  P = length(parameters(model))
  E = Int((P - 2N)/3)
  new_map_components = []
  for j in 1:N
    @nonamespace push!(new_map_components, (model.α[j], p_sample[i,j]))
    @nonamespace push!(new_map_components, (model.β[j], p_sample[i,N+j]./p_sample[i,N+1]))
  end
  for j in 1:E
    @nonamespace push!(new_map_components, (model.γ[j], p_sample[i,2*N+j]./p_sample[i,N+1]))
    @nonamespace push!(new_map_components, (model.κ[j], p_sample[i,2*N+E+j]))
    @nonamespace push!(new_map_components, (model.η[j], p_sample[i,2*N+2*E+j]))
  end
  new_pmap = Dict(new_map_components)
  new_p = ModelingToolkit.varmap_to_vars(new_pmap, parameters(model))
  return new_p
end

function generate_LHC_sample(model, samples, lim)
  N = length(species(model))
  P = length(parameters(model))
  E = Int((P - 2N)/3)
  LHCplan = randomLHC(samples, P)
  scaling_plan = Tuple{Float64,Float64}[]
  scaling_plan = vcat(scaling_plan, [(lim["α"][1], lim["α"][2]) for i in 1:N])
  scaling_plan = vcat(scaling_plan, [(lim["β"][1], lim["β"][2]) for i in 1:N])
  scaling_plan = vcat(scaling_plan, [(lim["γ"][1], lim["γ"][2]) for i in 1:E])
  scaling_plan = vcat(scaling_plan, [(lim["κ"][1], lim["κ"][2]) for i in 1:E])
  scaling_plan = vcat(scaling_plan, [(lim["η"][1], lim["η"][2]) for i in 1:E])
  p_sample = scaleLHC(LHCplan, scaling_plan)
  return p_sample
end

function create_generic_ode_problem(model)
  N = length(species(model))
  P = length(parameters(model))
  E = Int((P - 2N)/3)
  p_dict_components = []
  @nonamespace p_dict_components = vcat(p_dict_components, [(model.α[i], 0.5) for i in 1:N])
  @nonamespace p_dict_components = vcat(p_dict_components, [(model.β[i], 1.0) for i in 1:N])
  @nonamespace p_dict_components = vcat(p_dict_components, [(model.γ[i], 1.0) for i in 1:E])
  @nonamespace p_dict_components = vcat(p_dict_components, [(model.κ[i], 0.5) for i in 1:E])
  @nonamespace p_dict_components = vcat(p_dict_components, [(model.η[i], 1.0) for i in 1:E])
  pmap = Dict(p_dict_components)
  @nonamespace u0_dict_components = [(model.a[i], 0.5) for i in 1:N]
  u₀map = Dict(u0_dict_components)
  # Setup ODE problem
  prob = ODEProblem(model, u₀map, (0.0, 1.0), pmap)
  return prob
end

function infer_timescale(model, p_sample, i)
  N = length(species(model))
  P = length(parameters(model))
  E = Int((P - 2N)/3)
  rates = [p_sample[i,k] for k in (N+1):(2*N+E)]./p_sample[i,N+1]
  timescale = 1.0/minimum(rates)
  return timescale
end

function find_oscillations(model, samples, param_limits)
  # Periodogram hyperparameter
  f_sampling = 4000
  # Total simulation time = equil_tscales * inferred_tscale
  equil_tscales = 10
  # Create parameter sample with LHC
  N = length(species(model))
  P = length(parameters(model))
  E = Int((P - 2N)/3)
  p_sample = generate_LHC_sample(model, samples, param_limits)
  # Create generic parameter map and initial condition
  prob = create_generic_ode_problem(model)

  # Define functions to work with EnsembleProblem
  function prob_func(prob, i, repeat)
    i = Int(i)
    #=
    # Report status to terminal
    if 100*i/samples % 10 == 0
      println(100*i/samples)
    end
    =#
    # Set random initial condition
    equil_u0 = rand(N)
    # Select new parameters
    new_p = select_params_from_sample(model, p_sample, i)
    # Estimate timescale
    timescale = infer_timescale(model, p_sample, i)
    # Equilibrate
    equilibration = solve(prob, lsoda(), p=new_p, u0=equil_u0,
                          tspan=[0.0, equil_tscales*timescale])
    # Set new initial condition
    new_u0 = equilibration.u[end]
    remake(prob, p=new_p, u0=new_u0, tspan=[0.0, equil_tscales*timescale])
  end

  function output_func(sol, i)
    t_eval = LinRange(sol.t[1], sol.t[end], f_sampling)
    peaks = find_periodogram_peak(sol, t_eval, f_sampling)
    out = [sol, peaks]
    return (out, false)
  end
  # Setup ensemble problem and solve it
  ensemble_prob = EnsembleProblem(prob, prob_func=prob_func,
                                  output_func=output_func)
  sim = solve(ensemble_prob, lsoda(), EnsembleThreads(),
              trajectories=samples, abstol=1e-12)
  # Check if any solution oscillates
  osci_idxs = Int[]
  for i in 1:samples
    freq_comparison = [sim.u[i][2][j][1] == sim.u[i][2][j+1][1] for j in 1:N-1]
    if all(freq_comparison)
      peak_height = [sim.u[i][2][j][2] > 1e-6 for j in 1:N]
      if all(peak_height)
        push!(osci_idxs, i)
      end
    end
  end
  return [sim, p_sample, osci_idxs]
end
