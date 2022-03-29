function create_model(connectivity)
  N = size(connectivity, 1) # Number of nodes
  E = count(!iszero, connectivity)# Number of edges
  @parameters α[1:N], β[1:N], γ[1:E], κ[1:E], η[1:E]
  @variables t, a[1:N](t)
  rxs = Reaction[]
  for i=1:size(connectivity, 1)
    push!(rxs, Reaction(β[i]*α[i], nothing, [a[i]]))
    push!(rxs, Reaction(β[i], [a[i]], nothing))
  end
  e = 1 # Edge count
  for (i, row) in enumerate(eachrow(connectivity))
    for (j, val) in enumerate(row)
      if val == -1
        push!(rxs, Reaction(hill(abs(a[j]),γ[e],κ[e],η[e]), [a[i]], nothing))
        e += 1
      elseif val == 1
        push!(rxs, Reaction((1.0 - a[i]) * hill(abs(a[j]),γ[e],κ[e],η[e]), nothing, [a[i]]))
        e += 1
      end
    end
  end
  @named model = ReactionSystem(rxs, t)
  return model
end
