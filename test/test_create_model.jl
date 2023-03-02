"""
Evaluates if create_model returns a model with the correct number of nodes, parameters, and reactions.
Additionally, tests if proper errors are thrown when input is:
- An empty matrix 
- Not a two dimensional matrix
- Not a square matrix
- Has values that are not 0, 1, or -1
"""
# Connectivities resulting in functional models
connectivities = [[0 1; -1 0], [0 0; -1 0], [1 0; 0 0], 
                  [0 0 1; 1 0 0; 0 -1 0], [0 0 -1; -1 0 0; 0 -1 0], 
                  [0 0 1 0; 1 0 0 -1; 0 -1 0 0; 0 0 1 0], 
                  [0 0 0 -1; 0 -1 0 0; 0 0 -1 0; -1 1 0 1]] 
true_nodes = [2, 2, 2, 3, 3, 4, 4]
true_parameters = [10, 7, 7, 15, 15, 23, 26]
true_reactions = [8, 7, 7, 12, 12, 17, 18]
# Connectivities resulting in errors
error_connectivities = [
    [], [;;], zeros((3,3)), ones((1,2)), 
    [1 2 0; 0 0 0; 0 0 0], [1 0 0; 0 -2 0; 0 0 0]] 
expected_error_messages = [
    "Connectivity cannot be empty", "Connectivity cannot be empty",
    "Connectivity has to be a 2x2 matrix", "Connectivity has to be a square matrix",
    "Only -1, 0, and 1 are allowed as connectivity values",
    "Only -1, 0, and 1 are allowed as connectivity values"]

for (idx, c) in enumerate(connectivities)
    @test length(species(create_model(c))) == true_nodes[idx]
    @test length(parameters(create_model(c))) == true_parameters[idx]
    @test length(reactions(create_model(c))) == true_reactions[idx]
end

for (idx, c) in enumerate(error_connectivities)
    let err = nothing
        try
            create_model(c)
        catch err
            @test err isa DomainError
            @test occursin(expected_error_messages[idx], sprint(showerror, err))
        end
    end
end



