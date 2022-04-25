function peak_similarity(connectivity, pmap, tspan, init_cond)
  model = create_model(connectivity)
  @nonamespace u₀map = [model.a[1] => init_cond[1], model.a[2] => init_cond[2],
                        model.a[3] => init_cond[3]]
  prob = ODEProblem(model, u₀map, tspan, pmap)
  sol = solve(prob, lsoda())
  f_sampling = Int(2e3)
  t_eval = LinRange(tspan[1], tspan[end], f_sampling)
  fft_sampling = round(Int, f_sampling/tspan[end])
  pks = find_periodogram_peak(sol, t_eval, fft_sampling)
  return [pks[1][1] == pks[2][1], pks[2][1] == pks[3][1]]
end

# Goodwin Oscillator
connectivity = [0 0 1; 1 0 0; 0 -1 0]

# Test Goodwin's oscillatory solution
@nonamespace pmap_osci = (model.α[1] => 0.00004564, model.α[2] => 0.16342066,
                     model.α[3] => 0.70351496, model.β[1] => 56.67356362,
                     model.β[2] => 50.87885942, model.β[3] => 30.31904319,
                     model.γ[1] => 1657.39796806, model.γ[2] => 6988.78292809,
                     model.γ[3] => 7461.31128709, model.κ[1] => 0.11196140,
                     model.κ[2] => 0.50912113, model.κ[3] => 0.33210760,
                     model.η[1] => 3.67709117, model.η[2] => 7.17307667,
                     model.η[3] => 6.80130072)
init_cond_osci = [0.21, 0.23, 0.039]
# We expect all periodogram peaks to have the same frequency value for an
# oscillatory solution
@test all(peak_similarity(connectivity, pmap_osci, (0.0, 1.0), init_cond_osci))

# Test Goodwin's steady state solution
@nonamespace pmap_ss = (model.α[1] => 0.00004564, model.α[2] => 0.16342066,
                     model.α[3] => 0.70351496, model.β[1] => 6.67356362,
                     model.β[2] => 50.87885942, model.β[3] => 300.31904319,
                     model.γ[1] => 1657.39796806, model.γ[2] => 6988.78292809,
                     model.γ[3] => 7461.31128709, model.κ[1] => 0.11196140,
                     model.κ[2] => 0.50912113, model.κ[3] => 0.33210760,
                     model.η[1] => 3.67709117, model.η[2] => 7.17307667,
                     model.η[3] => 6.80130072)
init_cond_ss = [0.577, 0.99, 0.027]
# In a steady state solution we expect to see different frequencies for the
# peaks in the periodogram
@test !all(peak_similarity(connectivity, pmap_ss, (0.0, 1.0), init_cond_ss))
