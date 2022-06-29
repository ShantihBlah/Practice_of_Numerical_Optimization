include("problems/unconstraint/quadratic_func.jl")
# include("problems/unconstraint/Rosenbrock_func.jl")

include("Chapter_3/line_search.jl")
include("Chapter_4/trust_region.jl")


test_result = lineSearch(quadratic, test_vec, 100, BK_method="Newtown")
y_value = quadratic(test_result)

# test_result = trustRegion(quadratic, test_vec, 100, delta_max=20)
# y_value = quadratic(test_result)

# test_result = trustRegion(rosenbrock, test_vec, 100, delta_max=20)
# y_value = rosenbrock(test_result)

println("Initial X: ", test_vec)
println("Initial Y: ", init_y)
println("Final X: ", test_result)
println("Final Y: ", y_value)