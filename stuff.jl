using LinearAlgebra
using DifferentialEquations
using Plots

const WATER_RPM = 25
const WATER_RAD_VEL = WATER_RPM * 2π / 60 # 0.837758
# const M_PARTICLE = 1e-12 # kg
const M_PARTICLE = 1 * 70 / (3 * 10^13)
const R_PARTICLE = 5e-6
const A_PARTICLE = π * R_PARTICLE^2
const V_PARTICLE = (4/3) * π * R_PARTICLE^3 # m^3
const ρ_PARTICLE = M_PARTICLE / V_PARTICLE
const Cd_PARTICLE = 0.4 # assumed sphere
const ρ_FLUID = 1000 # kg / m^3
const GRAV_ACCEL = 9.81 # m/s^2
const R_WALLZ = 0.014 # 1.4 cm


# Pp is position vector of particle
function fluid_velocity(Pp::AbstractVector, ω=WATER_RAD_VEL)
    r = sqrt(sum(abs2, Pp))
    a_x, a_y = Pp[1], Pp[2]
    v_w_dir = normalize!([-a_y, a_x], 2) # [-a_y, a_x] for CCW, [a_y, -a_x] for CW
    v_w_magnitude = ω * r
    v_w = v_w_magnitude * v_w_dir
    return v_w
end

function drag_force_vec(p_particle, v_particle, ρ=ρ_FLUID, A=A_PARTICLE, Cd=Cd_PARTICLE)
    v_f = fluid_velocity(p_particle)
    v_rel_vec = v_f - v_particle
    v_rel = norm(v_rel_vec, 2)
    if iszero(v_rel)
        return [0.0, 0.0]
    else
        v_rel_dir = normalize!(v_rel_vec, 2)
        Fd = 0.5 * ρ * v_rel^2 * A * Cd
        Fd_vec = Fd * v_rel_dir
        return Fd_vec
    end
end

function buoyant_force_vec(V=V_PARTICLE ,ρ=ρ_FLUID)
    g = GRAV_ACCEL
    return [0.0, ρ * g * V]
end

function gravity_force_vec(m=M_PARTICLE)
    g = GRAV_ACCEL
    return [0.0, -m*g]
end

#=
Missing:
    Pressure Gradient force
    Shear induced lift force
=#
function resultant_force_vec(X_particle::AbstractVector)
    p_particle = X_particle[1:2]
    v_particle = X_particle[3:4]
    return drag_force_vec(p_particle, v_particle) .+ buoyant_force_vec() .+ gravity_force_vec()
end

# F = ma
# F = m v̇
# v̇ = F(x,v) / m

#=
    X = [px, py, vx, vy]
    for initial condition:
        vx, vy equals velocity of water
=#
function particle_dynamics(X::AbstractVector, p, t)
    F = resultant_force_vec(X)
    x = X[1:2]
    if norm(x,2) ≥ R_WALLZ
        ẋ = fluid_velocity(x)
        v̇ = [0.0, 0.0]
        return [ẋ; v̇]
    else
        v = X[3:4]
        ẋ = v
        v̇ = F ./ M_PARTICLE
        return [ẋ ; v̇]
    end
end

p_particle = [0.0, 0.012]
v_particle = fluid_velocity(p_particle)
# v_particle = [0.0, 0.0]

u0 = [p_particle; v_particle]
tspan = (0.0,10.0)
prob = ODEProblem(particle_dynamics, u0, tspan)
sol = solve(prob, Tsit5())

p = plot(getindex.(sol.u,1),getindex.(sol.u,2), xlabel="X (m)", ylabel="y (m)", label="", aspect_ratio=:equal)

# savefig(p, "try8rpm.png")

#
drag_force_vec(p_particle, v_particle)



drag_force_vec(p_particle, v_particle) / M_PARTICLE
buoyant_force_vec() / M_PARTICLE
gravity_force_vec() / M_PARTICLE


##
2 * R_WALLZ * (0.9*R_WALLZ) * WATER_RAD_VEL / (0.69 / ((10^3)^2)) 

using HARV
using StaticArrays
model = HARV.HARVModel()
dynamics = HARV.particle_dynamics(model)
u0 = [0.0, eps(), 0.0, 0.0]


using BenchmarkTools
dynamics(u0, nothing, nothing)
@benchmark dynamics(u0, nothing, nothing)

v1 = SA[1.,2.]
v2 = [2.,3.]
