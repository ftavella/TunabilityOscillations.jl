samples = Int(2e3)
"""
networks = Dict([
  ("G", [0 0 1; 1 0 0; 0 -1 0]),
  ("R", [0 0 -1; -1 0 0; 0 -1 0]),
  ("S", [1 -1; 1 0]),
  ("G+IPFC", [0 0 1; 1 0 0; 0 -1 1]),
  ("G+INFB", [0 0 1; 1 -1 0; 0 -1 0]),
  ("F", [-1 0 0 1; 1 -1 0 0; 0 1 -1 0; 0 0 -1 1]),
  ("G+R", [0 -1 1; 1 0 -1; -1 -1 0]),
  ("Cdc25", [0 0 -1 1; -1 0 0 0; 0 -1 0 0; 1 0 0 0]),
  ("Wee1", [0 0 1 1; 1 0 0 0; 0 -1 0 0; -1 0 0 0]),
])
"""
networks = Dict([
  ("G", [0 0 1; 1 0 0; 0 -1 0]),
  ("R", [0 0 -1; -1 0 0; 0 -1 0]),
])
hparams = Dict([("peak_num_thresh", 30), ("freq_tolerance", 0.01),
                ("equil_tscales", 50), ("sim_tscales", 40),
                ("power_threshold", 1e-7), ("amp_cv_thresh", 0.05),
                ("abstol", 1e-6), ("reltol", 1e-3),
                ("solver", solver), ("maxiters", 1e6),
                ("sampling_style", "lhc"), ("sampling_scale", "log")])
param_limits = Dict([
  ("α", [10^-2.5, 10^2.5]),
  ("β", [10^-2.5, 10^2.5]),
  ("γ", [10^2, 10^4]),
  ("κ", [0.1, 1.0]),
  ("η", [1.0, 10.0]),])

for (name, connectivity) in networks
  println("Analyzing network $name, Connectivity: $connectivity")
  local model = create_model(connectivity)
  @time sim, p_sample = find_oscillations(model, samples, param_limits, hparams)
  @test size(p_sample, 1) > 1
end

networks_not_oscil = Dict([
  ("P3", [0 0 1; 1 0 0; 0 1 0]),
  ("RP", [0 0 1; -1 0 0; 0 -1 0]),
  ("P2", [0 1; 1 0]),
  ("P4", [0 0 0 1; 1 0 0 0; 0 1 0 0; 0 0 1 0]),
  ("P5", [0 0 0 -1; 1 0 0 0; 0 -1 0 0; 0 0 1 0]),
])

for (name, connectivity) in networks_not_oscil
  println("Analyzing network $name, Connectivity: $connectivity")
  local model = create_model(connectivity)
  @time sim, p_sample = find_oscillations(model, samples, param_limits, hparams)
  @test size(p_sample, 1) == 0
end
