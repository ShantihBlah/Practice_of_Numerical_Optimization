using LinearAlgebra

function preliminaryConjugateGradient(A_matrix, x_vec, b_vec)
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
        # Reset init_point
        # if k % length(b_vec) == 0
        #     println(k)
        #     residual_k = A_matrix * x_vec_k - b_vec
        #     p_k = -residual_k
        # end
    end
    println("Iteration number: ", k)
    println("Residual: ", norm(A_matrix * x_vec_k - b_vec))
    return x_vec_k
end

function standardConjugateGradient(A_matrix, x_vec, b_vec)
    tolerance = 1e-16
    x_vec_k = x_vec
    residual_k = A_matrix * x_vec_k - b_vec
    p_k = -residual_k
    k = 0
    while norm(residual_k) > tolerance
        # println(k, ": ", norm(residual_k))
        # println(k, ": ", norm(A_matrix * x_vec_k - b_vec))
        residual_k_square = Transpose(residual_k)*residual_k
        quadratic_term = Transpose(p_k)*A_matrix*p_k
        alpha_k = residual_k_square / quadratic_term
        x_vec_k_1 = x_vec_k + alpha_k * p_k
        residual_k_1 = residual_k + alpha_k * A_matrix * p_k
        beta_k_1 = (Transpose(residual_k_1)*residual_k_1) / residual_k_square
        p_k_1 = -residual_k_1 + beta_k_1*p_k
        x_vec_k = x_vec_k_1
        p_k = p_k_1
        residual_k = residual_k_1
        k = k + 1
        # Reset init_point
        # if k % length(b_vec) == 0
        #     println(k)
        #     residual_k = A_matrix * x_vec_k - b_vec
        #     p_k = -residual_k
        # end
    end
    println("Iteration number: ", k)
    println("Residual: ", norm(A_matrix * x_vec_k - b_vec))
    return x_vec_k
end

