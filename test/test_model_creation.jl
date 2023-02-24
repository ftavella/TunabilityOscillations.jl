function test_model_creation(networks, simulated_peaks)
  """
  Evaluates if the function create_model returns a model with the correct number
  of species, parameters, and reactions.  Additionally, it simulates that model and 
  corroborates that solutions are periodic by matching the expected number of peaks. 
  The information for the reference networks is stored in  test/oscillatory_parameter_sets_*.jld2
  """
  for (name, connectivity) in networks
    println("Analyzing network $name, Connectivity: $connectivity")
    data = jldopen("oscillatory_parameter_sets_$name.jld2")
    connectivity = data["connectivity"]
    sim_result = data["sim_result"]
    p_sample = data["p_sample_result"]
    hparams = data["hparams"]
    model = create_model(connectivity)
    N = length(species(model))
    P = length(parameters(model))
    E = Int((P - 2N)/3)

    # Test that the model has the correct amount of species, parameters, and reactions
    @test length(species(model)) == size(connectivity, 1)
    @test length(parameters(model)) == 2*size(connectivity, 1) + 3*count(!iszero, connectivity)
    @test length(reactions(model)) == 3*size(connectivity, 1) + count(!iszero, connectivity)

    for (idx, p_set) in enumerate(eachrow(p_sample))
        period = 1.0 / sim_result[idx][3][1]
        initial_state = sim_result[idx][2]

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
        ode_params = ModelingToolkit.varmap_to_vars(Dict(new_map_vals), parameters(model))

        p_dict_vals = []
        @nonamespace p_dict_vals = vcat(p_dict_vals, [(model.α[i], 0.5) for i in 1:N])
        @nonamespace p_dict_vals = vcat(p_dict_vals, [(model.β[i], 1.0) for i in 1:N])
        @nonamespace p_dict_vals = vcat(p_dict_vals, [(model.γ[i], 1.0) for i in 1:E])
        @nonamespace p_dict_vals = vcat(p_dict_vals, [(model.κ[i], 0.5) for i in 1:E])
        @nonamespace p_dict_vals = vcat(p_dict_vals, [(model.η[i], 1.0) for i in 1:E])
        @nonamespace u0_dict_vals = [(model.a[i], 0.5) for i in 1:N]
        pmap = Dict(p_dict_vals)
        u₀map = Dict(u0_dict_vals)
        prob = ODEProblem(model, u₀map, (0.0, 1.0), pmap)
        prob = remake(prob, p=ode_params, u0=initial_state, tspan=[0.0, simulated_peaks*period])
        sol = solve(prob, hparams["solver"], abstol=1e-7, reltol=1e-4, maxiters=hparams["maxiters"], 
                    dense=false)
        
        # Is solution is oscillatory?
        rough_amp = maximum(sol[1,:]) - minimum(sol[1,:])
        pks_max, ~ = findmaxima(sol[1,:])
        pks_max, ~ = peakproms(pks_max, sol[1,:], minprom=rough_amp/3.0)

        @test size(pks_max, 1) >= simulated_peaks - 2

    end
  end
end

# Networks to be tested
networks = Dict([
  ("G", [0 0 1; 1 0 0; 0 -1 0]),
  ("G+NAA+NBB+PCC", [-1 0 1; 1 -1 0; 0 -1 1]),
  ("R", [0 0 -1; -1 0 0; 0 -1 0]),
  ("R+PAA+PBB+PCC", [1 0 -1; -1 1 0; 0 -1 1]),
])
# Number of peaks to simulate
simulated_peaks = 5
test_model_creation(networks, simulated_peaks)
