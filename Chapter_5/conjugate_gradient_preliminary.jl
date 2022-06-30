using LinearAlgebra

function conjugate_gradient_preliminary(A_matrix, x_vec, b_vec)
    tolerance = 1e-16
    x_vec_k = x_vec
    residual_k = A_matrix * x_vec_k - b_vec
    p_k = -residual_k
    k = 0
    while norm(residual_k) > tolerance
        # println(k, ": ", norm(residual_k))
        quadratic_term = (Transpose(p_k)*A_matrix*p_k)
        alpha_k = - (Transpose(residual_k)*p_k) / quadratic_term
        x_vec_k_1 = x_vec_k + alpha_k * p_k
        residual_k_1 = A_matrix * x_vec_k_1 - b_vec
        beta_k_1 = (Transpose(residual_k_1)*A_matrix*p_k) / quadratic_term
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
