#
#  Solve A*x = b
#

using Random

linear_transform_rng = MersenneTwister(1234)
linear_transform_dim = 100
linear_transform_B = randexp(linear_transform_rng, linear_transform_dim, linear_transform_dim)
linear_transform_x = 5.0*ones(linear_transform_dim)
linear_transform_b = zeros(linear_transform_dim)
linear_transform_A = Transpose(linear_transform_B) * linear_transform_B

linear_transform_setup = Dict("A" => linear_transform_A, "x" => linear_transform_x, "b" => linear_transform_b)