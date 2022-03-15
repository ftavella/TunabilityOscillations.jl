function create_model_expression(connectivity)
  # Create and aggregate interaction terms based on the connectivity
  interaction_terms = Expr[]
  for (i, row) in enumerate(eachrow(connectivity))
    row_terms = String[]
    for (j, val) in enumerate(row)
      if val != 0
        # Define the substrate and sign of the interaction
        sub = (val == 1) ? "+ (1.0 - abs(x[$i]))" : "- abs(x[$i])"
        numerator = "$sub * γ[$i, $j] * abs(x[$j]) ^ η[$i, $j]"
        denominator = "(abs(x[$j]) ^ η[$i, $j] + κ[$i, $j] ^ η[$i, $j]) "
        term = numerator * " / " * denominator
        push!(row_terms, term)
      end
    end
    push!(interaction_terms, Meta.parse(join(row_terms)))
  end
  # Create full terms for each variable
  expr_terms = Expr[]
  for i in 1:size(connectivity, 1)
    t = :(dx[$i] = β[$i] * (α[$i] - abs(x[$i])) +  $(interaction_terms[i]))
    push!(expr_terms, t)
  end
  # Define the expression of the model's equations
  expr = quote
    function model_equations!(dx, x, p, t)
      (α, β, γ, κ, η) = p
      $(expr_terms...)
    end
  end

  return expr
end
