samples = Int(1e3)
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
param_limits = Dict([
  ("α", [0.0, 1.0]),
  ("β", [10.0, 100.0]),
  ("γ", [100.0, 10000.0]),
  ("κ", [0.0, 1.0]),
  ("η", [3.0, 10.0]),
])
for (name, connectivity) in networks
  println("Analyzing network $name, connectivity: $connectivity")
  local model = create_model(connectivity)
  sim, p_sample, osci_idxs = find_oscillations(model, samples, param_limits)
  @test length(osci_idxs) > 1
end

networks_not_oscil = Dict([
  ("P3", [0 0 1; 1 0 0; 0 1 0]),
  ("RP", [0 0 1; -1 0 0; 0 -1 0]),
  ("N2", [0 -1; 1 0]),
  ("P4", [0 0 0 1; 1 0 0 0; 0 1 0 0; 0 0 1 0]),
  ("P5", [0 0 0 -1; 1 0 0 0; 0 -1 0 0; 0 0 1 0]),
])

for (name, connectivity) in networks_not_oscil
  println("Analyzing network $name, connectivity: $connectivity")
  local model = create_model(connectivity)
  sim, p_sample, osci_idxs = find_oscillations(model, samples, param_limits)
  @test length(osci_idxs) == 0
end

  # scaling_plan = vcat(scaling_plan, [(0.0, 1.0) for i in 1:N]) # α
  # scaling_plan = vcat(scaling_plan, [(10.0, 100.0) for i in 1:N]) # β
  # scaling_plan = vcat(scaling_plan, [(100.0, 10000.0) for i in 1:E]) # γ
  # scaling_plan = vcat(scaling_plan, [(0.0, 1.0) for i in 1:E]) # κ
  # scaling_plan = vcat(scaling_plan, [(3.0, 10.0) for i in 1:E]) # η
