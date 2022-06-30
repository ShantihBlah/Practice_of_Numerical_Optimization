using Random

include("problems/unconstraint/quadratic_func.jl")
include("problems/unconstraint/Rosenbrock_func.jl")

include("Chapter_3/line_search.jl")
include("Chapter_4/trust_region.jl")
include("Chapter_5/conjugate_gradient_preliminary.jl")
include("Chapter_5/conjugate_gradient.jl")

# problem_setup = quadratic_setup
# # problem_setup = rosenbrock_setup

# object_function = problem_setup["obj_func"]
# init_vec = problem_setup["init_vec"]
# init_y = object_function(init_vec)
# max_iter_num = 2000

# @time solution_vec = lineSearch(object_function, init_vec, max_iter_num, BK_method="Speedest_Descend")
# # @time solution_vec = trustRegion(object_function, init_vec, max_iter_num, delta_max=20)

# solution_y = object_function(solution_vec)

# println("Initial X: ", init_vec)
# println("Initial Y: ", init_y)
# println("Final X: ", solution_vec)
# println("Final Y: ", solution_y)


rng = MersenneTwister(1234)

dim = 100
B = randexp(rng, dim, dim)
x = 5.0*ones(dim)
b = zeros(dim)
A = (Transpose(B) * B)


@time solution_x = conjugate_gradient_preliminary(A, x, b)
@time solution_x = conjugate_gradient(A, x, b)
