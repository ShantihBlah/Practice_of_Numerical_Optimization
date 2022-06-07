# include("problems/unconstraint/quadratic_func.jl")
include("problems/unconstraint/Rosenbrock_func.jl")

include("Chapter_3/line_search.jl")
include("Chapter_4/trust_region.jl")


# test_result = lineSearch(test_obj_func, test_diff_func, test_Hessian, test_vec, 100, BK_method="Newtown")
test_result = trustRegion(test_obj_func, test_diff_func, test_Hessian, test_vec, 1000, delta_max=1)

y_value = test_obj_func(test_result)

println("Initial X: ", test_vec)
println("Initial Y: ", init_y)
println("Final X: ", test_result)
println("Final Y: ", y_value)