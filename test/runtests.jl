using Test
using DifferentialEquations
using HARV

@testset "smoke" begin
    model = HARVModel()
    dynamics = particle_dynamics(model)
    u0 = [0.0, eps(), 0.0, 0.0]
    tspan = (0., 1.)
    prob = ODEProblem(dynamics, u0, tspan)
    sol = solve(prob, Tsit5())
    @test sol isa ODESolution
end
