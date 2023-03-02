"""
    create_model(connectivity::Matrix)

Creates a ReactionSystem based on the provided `connectivity`.

# Arguments
- `connectivity::Matrix`: A 2 dimensional matrix filled with -1, 0, and 1 values indicating the edges of the network.
"""
function create_model(connectivity::Matrix)
  # Check that input is correct
  if isempty(connectivity)
    throw(DomainError(connectivity, "Connectivity cannot be empty"))
  elseif ndims(connectivity) != 2
    throw(DomainError(connectivity, "Connectivity has to be a 2x2 matrix"))
  elseif size(connectivity, 1) != size(connectivity, 2)
    throw(DomainError(connectivity, "Connectivity has to be a square matrix"))
  elseif any([v ∉ [-1 0 1] for v in connectivity])
    throw(DomainError(connectivity, "Only -1, 0, and 1 are allowed as connectivity values"))
  end
  # Number of nodes
  N = size(connectivity, 1) 
  # Number of edges
  E = count(!iszero, connectivity)
  @parameters α[1:N], β[1:N], γ[1:E], κ[1:E], η[1:E]
  @variables t, a[1:N](t)
  rxs = Reaction[]
  for i=1:size(connectivity, 1)
    # + Alpha
    push!(rxs, Reaction(α[i], nothing, [a[i]]))
    # - Alpha x A
    push!(rxs, Reaction(α[i], [a[i]], nothing))
    # - Beta x A
    push!(rxs, Reaction(β[i], [a[i]], nothing))
  end
  # Edge count
  e = 1 
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
