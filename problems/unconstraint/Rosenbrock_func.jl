using LinearAlgebra


function test_obj_func(vec::Vector)
    return 100*(vec[2]-vec[1]^2)^2 + (1-vec[1])^2
end


function test_diff_func(vec::Vector)
    diff_vec = [400*vec[1]*(vec[1]^2-vec[2]) + 2*(vec[1]-1); 200*(vec[2]-vec[1]^2)]
    return diff_vec
end

function test_Hessian(vec::Vector)
    Q = [400*(3*vec[1]^2-vec[2])+2 -400*vec[1]; -400*vec[1] 200]
    return Q
end

test_vec = [4.;2.]
init_y = test_obj_func(test_vec)