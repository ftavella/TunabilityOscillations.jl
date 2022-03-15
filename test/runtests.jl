using TunabilityOscillations
using DifferentialEquations
using Test

using Plots

@testset "TunabilityOscillations.jl" begin
    # Test that the created model equations behave like models defined by hand
    function Goodwin!(du,u,p,t)
      x, y, z = u
      (α, β, γ, κ, η) = p
      du[1] = β[1]*(α[1] - x) + γ[1,3] * (1.0 - x) * abs(z)^η[1,3] / (abs(z)^η[1,3] + κ[1,3]^η[1,3])
      du[2] = β[2]*(α[2] - y) + γ[2,1] * (1.0 - y) * abs(x)^η[2,1] / (abs(x)^η[2,1] + κ[2,1]^η[2,1])
      du[3] = β[3]*(α[3] - z) - γ[3,2] * z * abs(y)^η[3,2] / (abs(y)^η[3,2] + κ[3,2]^η[3,2])
    end
    connectivity_goodwin = [0 0 1; 1 0 0; 0 -1 0]
    expr_goodwin = create_model_expression(connectivity_goodwin)
    eval(expr_goodwin)

    # Parameters Goodwin
    u0 = [0, 0, 0]
    α = [0.00004564, 0.16342066, 0.70351496]
    β = [56.67356362, 50.87885942, 30.31904319]
    γ = [0 0 1657.39796806
        6988.78292809 0 0
        0 7461.31128709 0]
    κ = [0 0 0.11196140
        0.50912113 0 0
        0 0.33210760 0]
    η = [0 0 3.67709117
        7.17307667 0 0
        0 6.80130072 0]
    p = [α, β, γ, κ, η]

    u0_values = rand(1000, 3)
    result_expr = zeros(1000, 3)
    result_fun = zeros(1000, 3)
    for (i, row) in enumerate(eachrow(u0_values))
      du_fun = zeros(1, 3)
      du_expr = zeros(1, 3)
      Goodwin!(du_fun, row, p, 0.0)
      model_equations!(du_expr, row, p, 0.0)
      result_fun[i,:] = du_fun
      result_expr[i,:] = du_expr
    end

    @test isapprox(result_expr, result_fun)

    function Goodwin_PF!(du, u,p,t)
      x, y, z = u
      (α, β, γ, κ, η) = p
      du[1] = β[1]*(α[1] - x) + γ[1,3] * (1.0 - x) * abs(z)^η[1,3] / (abs(z)^η[1,3] + κ[1,3]^η[1,3])
      du[2] = β[2]*(α[2] - y) + γ[2,1] * (1.0 - y) * abs(x)^η[2,1] / (abs(x)^η[2,1] + κ[2,1]^η[2,1])
      du[3] = β[3]*(α[3] - z) - γ[3,2] * z * abs(y)^η[3,2] / (abs(y)^η[3,2] + κ[3,2]^η[3,2]) + γ[3,3] * (1.0 - z) * abs(z)^η[3,3] / (abs(z)^η[3,3] + κ[3,3]^η[3,3])
    end

    connectivity_goodwin_pf = [0 0 1; 1 0 0; 0 -1 1]
    expr_goodwin_pf = create_model_expression(connectivity_goodwin_pf)
    eval(expr_goodwin_pf)

    # Parameters Goodwin + PF
    u0 = [0.5, 0.5, 0.5]
    α = [0.00012570, 0.24512445, 0.74518621]
    β = [70.08309071, 30.10686277, 37.80616599]
    γ = [0 0 6451.62448135
        1104.67978494 0 0
        0 7529.60087920 5429.62993215]
    κ = [0 0 0.66066269
        0.59212769 0 0
        0 0.69523907 0.58345202]
    η = [0 0 9.22476250
        4.34775355 0 0
        0 8.88820688 5.73136565]
    p = [α, β, γ, κ, η]

    u0_values = rand(1000, 3)
    result_expr = zeros(1000, 3)
    result_fun = zeros(1000, 3)
    for (i, row) in enumerate(eachrow(u0_values))
      du_fun = zeros(1, 3)
      du_expr = zeros(1, 3)
      Goodwin_PF!(du_fun, row, p, 0.0)
      model_equations!(du_expr, row, p, 0.0)
      result_fun[i,:] = du_fun
      result_expr[i,:] = du_expr
    end

    @test isapprox(result_expr, result_fun)

end
