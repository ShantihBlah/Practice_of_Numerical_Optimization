using LinearAlgebra
using ForwardDiff: gradient, hessian

function backtrackingLineSearch(object_function::Function, x_vec::Vector, p_vec::Vector)
    alpha_init = 1.0
    c = 0.001
    rho = 0.75
    alpha = alpha_init
    x_vec_k1 = x_vec .+ alpha.*p_vec
    # println(x_vec)
    while (object_function(x_vec_k1) > (object_function(x_vec) + c*alpha*Transpose(gradient(object_function, x_vec))*p_vec))
        x_vec_k1 = x_vec .+ alpha.*p_vec
        alpha *= rho
    end
    return alpha
end

