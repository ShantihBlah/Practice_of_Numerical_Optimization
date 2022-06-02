include("line_search.jl")

function test_obj_func(vec::Vector)
    A = [1. 0. 2.; 1. 2. 0.; 5. 0. 2.]
    Q = A * Transpose(A)
    return 0.5*transpose(vec)*Q*vec
end


function test_diff_func(vec::Vector)
    A = [1. 0. 2.; 1. 2. 0.; 5. 0. 2.]
    Q = A * Transpose(A)
    return Q*vec
end

function test_Hessian(vec::Vector)
    A = [1. 0. 2.; 1. 2. 0.; 5. 0. 2.]
    Q = A * Transpose(A)
    return Q
end

test_vec = [10.;4.;2.]
init_y = test_obj_func(test_vec)

# function test_obj_func(vec::Vector)
#     return 100*(vec[2]-vec[1]^2)^2 + (1-vec[1])^2
# end


# function test_diff_func(vec::Vector)
#     diff_vec = [400*vec[1]*(vec[1]^2-vec[2]) + 2*(vec[1]-1); 200*(vec[2]-vec[1]^2)]
#     return diff_vec
# end

# function test_Hessian(vec::Vector)
#     Q = [400*(3*vec[1]^2-vec[2])+2 -400*vec[1]; -400*vec[1] 200]
#     return Q
# end

# test_vec = [4.;2.]
# init_y = test_obj_func(test_vec)



# println(backtrackingLineSearch(test_obj_func, test_diff_func, test_vec, -test_diff_func(test_vec)))

# test_result = lineSearch(test_obj_func, test_diff_func, test_Hessian, test_vec, 40000, BK_method="Newtown")
test_result = lineSearch(test_obj_func, test_diff_func, test_Hessian, test_vec, 100, BK_method="Newtown")

y_value = test_obj_func(test_result)

println("Initial X: ", test_vec)
println("Initial Y: ", init_y)
println("Final X: ", test_result)
println("Final Y: ", y_value)