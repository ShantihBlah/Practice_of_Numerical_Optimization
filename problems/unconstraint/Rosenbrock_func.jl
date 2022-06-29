#
#  Rosenbrock function: https://en.wikipedia.org/wiki/Rosenbrock_function
#

using LinearAlgebra

function rosenbrock(x)
    a = one(eltype(x))
    b = 100 * a
    result = zero(eltype(x))
    for i in 1:length(x)-1
        result += (a - x[i])^2 + b*(x[i+1] - x[i]^2)^2
    end
    return result
end

test_vec = [4.;2.]
init_y = rosenbrock(test_vec)


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