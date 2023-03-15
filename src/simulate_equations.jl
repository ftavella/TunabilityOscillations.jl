"""
    reference_hparams::Dict

Reference hyperparameters used throughout the package

```@eval
reference_hparams = Dict("amp_cv_thresh" => 0.05, "sampling_style" => "lhc", "maxiters" => 1.0e7, 
                         "solver" => RadauIIA5(), "sampling_scale" => "log", "equil_tscales" => 10, 
                         "sim_tscales" => 20, "power_threshold" => 1.0e-7, "peak_num_thresh" => 15,
                         "abstol" => 1.0e-7, "freq_tolerance"  => 0.01, "reltol" => 0.0001)
display(reference_hparams)
```
"""
reference_hparams = Dict("amp_cv_thresh" => 0.05, "sampling_style" => "lhc", "maxiters" => 1.0e7, 
                         "solver" => RadauIIA5(), "sampling_scale" => "log", "equil_tscales" => 10, 
                         "sim_tscales" => 20, "power_threshold" => 1.0e-7, "peak_num_thresh" => 15,
                         "abstol" => 1.0e-7, "freq_tolerance"  => 0.01, "reltol" => 0.0001)

"""
    ode_params_from_param_set(model::ReactionSystem, p_set::AbstractVector)

Convert a vector of parameters to an array suitable for `ODEProblem` remaking using the `model`'s information.

# Arguments
- `model::ReactionSystem`: Model to be used as a reference.
- `p_set::AbstractVector`: Parameter values ordered as [α1, α2, ..., β1, ..., γ1, ..., κ1, ..., η1, ..., ηN] where N is the number of nodes.
"""
function ode_params_from_vector(model::ReactionSystem, p_set::AbstractVector)
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
    generic_ode_problem(model::ReactionSystem)

Create a generic ODE problem from a `model`.

# Arguments
- `model::ReactionSystem`: Model to be used as a reference.
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
    simulate(model::ReactionSystem, p_set::AbstractVector, integration_time::Float64, initial_condition::AbstractVector; hyperparameters::Dict=reference_hparams)

Integrate the differential equations specified by the `model` under the parameters on `p_set`. Returns the calculated solution. 

# Arguments
- `model::ReactionSystem`: Model to be used as a reference.
- `p_set::AbstractVector`: Parameter values ordered as [α1, α2, ..., β1, ..., γ1, ..., κ1, ..., η1, ..., ηN] where N is the number of nodes.
- `integration_time::Float64`: Total integration time
- `initial_condition::AbstractVector`: Initial state of each variable in the `model`
- `hyperparameters::Dict`: Solver hyperparameters. Defaults to reference parameters specified in `reference_hparams`
"""
function simulate(model::ReactionSystem, p_set::AbstractVector, integration_time::Float64, initial_condition::AbstractVector; hyperparameters::Dict=reference_hparams)
    ode_params = ode_params_from_vector(model, p_set)
    prob = generic_ode_problem(model)
    prob = remake(prob, p=ode_params, u0=initial_condition, tspan=[0.0, integration_time])
    sol = solve(prob, hyperparameters["solver"], abstol=hyperparameters["abstol"], 
                reltol=hyperparameters["reltol"], maxiters=hyperparameters["maxiters"])
    return sol
end