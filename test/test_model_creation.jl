function is_correct_model(connectivity, pmap, expected_final_state, tspan, init_cond)
  #=
    Function to test whether a generated model has the right number of species,
    parameters, reactions. This function also tests if the model reaches a
    desired final state when simulated from the provided initial conditions and
    time span.
  =#
  N = size(connectivity, 1) # Number of nodes
  E = count(!iszero, connectivity)# Number of edges

  model = create_model(connectivity)

  @test length(species(model)) == N
  @test length(parameters(model)) == 2*N + 3*E
  @test length(reactions(model)) == 2*N + E

  @nonamespace u₀map = [model.a[1] => init_cond[1], model.a[2] => init_cond[2],
                        model.a[3] => init_cond[3]]
  prob = ODEProblem(model, u₀map, tspan, pmap)
  sol = solve(prob, lsoda())
  error = broadcast(abs, sol.u[end,:][1] - expected_final_state)
  @test all(error .< 1e-2)
end

# Test Goodwin's oscillatory network
connectivity_g = [0 0 1; 1 0 0; 0 -1 0]
tspan_g = (0.0, 1.0)
init_cond_g = [0.9, 0.9, 0.9]
expected_final_state_g = [0.21, 0.23, 0.039]
model = create_model(connectivity_g)
@nonamespace pmap_g = (model.α[1] => 0.00004564, model.α[2] => 0.16342066,
                       model.α[3] => 0.70351496, model.β[1] => 56.67356362,
                       model.β[2] => 50.87885942, model.β[3] => 30.31904319,
                       model.γ[1] => 1657.39796806, model.γ[2] => 6988.78292809,
                       model.γ[3] => 7461.31128709, model.κ[1] => 0.11196140,
                       model.κ[2] => 0.50912113, model.κ[3] => 0.33210760,
                       model.η[1] => 3.67709117, model.η[2] => 7.17307667,
                       model.η[3] => 6.80130072)
is_correct_model(connectivity_g, pmap_g, expected_final_state_g, tspan_g, init_cond_g)

# Test Goodwin + PF
connect_gpf = [0 0 1; 1 0 0; 0 -1 1]
tspan_gpf = (0.0, 1.0)
init_cond_gpf = [0.9, 0.9, 0.9]
expec_final_state_gpf = [0.047, 0.43, 0.13]
model = create_model(connect_gpf)
@nonamespace pmap_gpf = (model.α[1] => 0.00012570, model.α[2] => 0.24512445,
                         model.α[3] => 0.74518621, model.β[1] => 70.08309071,
                         model.β[2] => 30.10686277, model.β[3] => 37.80616599,
                         model.γ[1] => 6451.62448135, model.γ[2] => 104.67978494,
                         model.γ[3] => 7529.60087920, model.γ[4] => 5429.62993215,
                         model.κ[1] => 0.66066269, model.κ[2] => 0.59212769,
                         model.κ[3] => 0.69523907, model.κ[4] => 0.58345202,
                         model.η[1] => 9.22476250, model.η[2] => 4.34775355,
                         model.η[3] => 8.88820688, model.η[4] => 5.73136565)
is_correct_model(connect_gpf, pmap_gpf, expec_final_state_gpf, tspan_gpf, init_cond_gpf)
