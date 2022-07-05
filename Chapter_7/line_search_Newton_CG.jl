using LinearAlgebra
using ForwardDiff: gradient, hessian

include("../Chapter_3/backtracking_line_search.jl")

function lineSearchNewtonCG(object_function::Function, x_vec::Vector, max_iter_num::Int, max_inner_iter_num::Int)
    tolerance = 1e-16
    x_vec_k = x_vec

    # setup the 1st Gradient
    f = object_function
    g = x_vec -> gradient(object_function, x_vec)
    H = x_vec -> hessian(object_function, x_vec)

    for i in range(0, length=max_iter_num)
        
        g_k = g(x_vec_k)
        B_k = H(x_vec_k)
        tolerance_k = min(0.5, sqrt(norm(g_k))) * norm(norm(g_k))
        z_j = zeros(length(x_vec_k))
        r_j = g(x_vec_k)
        d_j = -r_j

        p_k = zeros(length(x_vec_k))
        for j in range(0, length=max_inner_iter_num)
            sec_term = Transpose(d_j) * B_k * d_j
            if sec_term <= 0
                if j == 0
                    p_k = -g_k
                    break
                else
                    p_k = z_j
                    break
                end
            end
            alpha_j = (Transpose(r_j)*r_j) / (Transpose(d_j)*B_k*d_j)
            z_j = z_j + alpha_j*d_j
            beta_denominator = Transpose(r_j)*r_j
            r_j = r_j + alpha_j*B_k*d_j
            p_k = z_j
            if norm(r_j) < tolerance_k
                break
            end
            beta_j = (Transpose(r_j)*r_j) / beta_denominator
            d_j = -r_j + beta_j*d_j
        end
        p_vec = p_k
        alpha = backtrackingLineSearch(f, x_vec_k, p_vec)
        x_vec_k = x_vec_k + alpha * p_vec
        diff_value = norm(g(x_vec))
        if diff_value < tolerance
            break
        end
    end
    return x_vec_k
end