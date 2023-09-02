module HARV

using StaticArrays
using LinearAlgebra
using DifferentialEquations

include("model.jl")
export HARVModel, particle_dynamics

end # module HARV
