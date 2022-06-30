using LinearAlgebra
using Random

dim = 50

function quadratic(vec::Vector)
    rng = MersenneTwister(1234)
    A = randexp(rng, dim, dim)
    # A = [1. 0. 2.; 1. 2. 0.; 5. 0. 2.]
    Q = A * Transpose(A)
    result = 0.5*transpose(vec)*Q*vec
    return result
end

init_vec = 5.0*ones(dim)
# init_vec = [10.;4.;2.]

quadratic_setup = Dict("obj_func" => quadratic, "init_vec" => init_vec)