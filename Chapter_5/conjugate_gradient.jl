using LinearAlgebra
using ForwardDiff: gradient, hessian

function conjugate_gradient(A_matrix, x_vec, b_vec)
    tolerance = 1e-16
    x_vec_k = x_vec
    residual_k = A_matrix * x_vec_k - b_vec
    p_k = -residual_k
    k = 0
    while norm(residual_k) > tolerance
        # println(k, ": ", norm(residual_k))
        alpha_k = (Transpose(residual_k)*residual_k) / (Transpose(p_k)*A_matrix*p_k)
        x_vec_k_1 = x_vec_k + alpha_k * p_k
        residual_k_1 = residual_k + alpha_k * A_matrix * p_k
        beta_k_1 = (Transpose(residual_k_1)*residual_k_1) / (Transpose(residual_k)*residual_k) 
        p_k_1 = -residual_k_1 + beta_k_1*p_k
        x_vec_k = x_vec_k_1
        p_k = p_k_1
        residual_k = residual_k_1
        k = k + 1
        # if k == length(b_vec)
        #     break
        # end
    end
    println(k)
    return x_vec_k
end

