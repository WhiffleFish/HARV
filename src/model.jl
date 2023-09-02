const GRAV_ACCEL = 9.81

Base.@kwdef struct HARVModel
    ω_water::Float64 = 25 * 2π / 60 # 25 RPM by default
    m_particle::Float64 = 1 * 70 / (3 * 10^13)
    r_particle::Float64 = 5e-6
    a_particle::Float64 = π * r_particle^2
    v_particle::Float64 = (4/3) * π * r_particle^3
    Cd_particle::Float64 = 0.4 # assumed sphere
    ρ_fluid::Float64 = 1000 # kg/m^3
    r_walls::Float64 = 0.014 # 1.4 cm
end

fluid_velocity(model::HARVModel, p_particle) = fluid_velocity(p_particle, model.ω_water)

function fluid_velocity(p_particle::AbstractVector, ω)
    r = norm(p_particle, 2)
    a_x, a_y = p_particle[1], p_particle[2]
    v_w_dir = normalize(SA[-a_y, a_x], 2) # [-a_y, a_x] for CCW, [a_y, -a_x] for CW
    v_w_magnitude = ω * r
    v_w = v_w_magnitude * v_w_dir
    return v_w
end

function drag_force_vec(model::HARVModel, p_particle, v_particle)
    ρ = model.ρ_fluid
    A = model.a_particle
    Cd = model.Cd_particle

    v_f = fluid_velocity(model, p_particle)
    v_rel_vec = v_f - v_particle
    v_rel = norm(v_rel_vec, 2)
    if iszero(v_rel)
        return SA[0.0, 0.0]
    else
        v_rel_dir = normalize(v_rel_vec, 2)
        Fd = 0.5 * ρ * v_rel^2 * A * Cd
        Fd_vec = Fd * v_rel_dir
        return Fd_vec
    end
end

buoyant_force_vec(model::HARVModel) = buoyant_force_vec(model.v_particle, model.ρ_fluid)

function buoyant_force_vec(V::Float64, ρ::Float64)
    g = GRAV_ACCEL
    return SA[0.0, ρ * g * V]
end

gravity_force_vec(model::HARVModel) = gravity_force_vec(model.m_particle)

function gravity_force_vec(m::Float64)
    g = GRAV_ACCEL
    return SA[0.0, -m*g]
end

function resultant_force_vec(model::HARVModel, X_particle::AbstractVector)
    p_particle = X_particle[1:2]
    v_particle = X_particle[3:4]
    return drag_force_vec(model, p_particle, v_particle) .+ buoyant_force_vec(model) .+ gravity_force_vec(model)
end

#=
    X = [px, py, vx, vy]
    for initial condition:
        vx, vy equals velocity of water
=#

particle_dynamics(model::HARVModel) = (u,p,t) -> particle_dynamics(model, u, p, t)

function particle_dynamics(model::HARVModel, X::AbstractVector, p, t)
    F = resultant_force_vec(model, X)
    x = X[1:2]
    if norm(x,2) ≥ model.r_walls
        ẋ = fluid_velocity(model, x)
        v̇ = [0.0, 0.0]
        return [ẋ; v̇]
    else
        v = X[3:4]
        ẋ = v
        v̇ = F ./ model.m_particle
        return [ẋ ; v̇]
    end
end
