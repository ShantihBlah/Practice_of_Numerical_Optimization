#
#  min y(x) = x^T*A*x
#

using LinearAlgebra
using Random

quadratic_dim = 5

function quadratic(vec::Vector)
    rng = MersenneTwister(1234)
    A = randexp(rng, quadratic_dim, quadratic_dim)
    # A = [1. 0. 2.; 1. 2. 0.; 5. 0. 2.]
    Q = A * Transpose(A)
    result = 0.5*transpose(vec)*Q*vec
    return result
end

quadratic_init_vec = 5.0*ones(quadratic_dim)
# init_vec = [10.;4.;2.]

quadratic_setup = Dict("obj_func" => quadratic, "init_vec" => quadratic_init_vec)