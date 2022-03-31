samples = Int(1e3)
networks = Dict([
  ("G", [0 0 1; 1 0 0; 0 -1 0]),
  ("R", [0 0 -1; -1 0 0; 0 -1 0]),
  ("G+IPFC", [0 0 1; 1 0 0; 0 -1 1]),
  ("G+INFB", [0 0 1; 1 -1 0; 0 -1 0]),
  ("F", [-1 0 0 1; 1 -1 0 0; 0 1 -1 0; 0 0 -1 1]),
  ("G+R", [0 -1 1; 1 0 -1; -1 -1 0]),
  ("Cdc25", [0 0 -1 1; -1 0 0 0; 0 -1 0 0; 1 0 0 0]),
  ("Wee1", [0 0 1 1; 1 0 0 0; 0 -1 0 0; -1 0 0 0]),
])
connectivity = networks["G"]
model = create_model(connectivity)
sim, osci_idxs = find_oscillations(model, samples)
@test length(osci_idxs) > 1
