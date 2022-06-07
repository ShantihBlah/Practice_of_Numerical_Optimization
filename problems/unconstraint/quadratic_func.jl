using LinearAlgebra

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