using LinearAlgebra
using ForwardDiff: gradient, hessian

include("../Chapter_3/backtracking_line_search.jl")

function lineSearchBFGS(object_function::Function, x_vec::Vector, max_iter_num::Int)
    tolerance = 1e-16
    x_vec_current = x_vec
    H_0 = Diagonal(ones(size(x_vec)[1]))
    H_k = H_0
    x_vec_previous = x_vec_current

    # setup the 1st Gradient
    g = x_vec -> gradient(object_function, x_vec)

    for i in range(1, length=max_iter_num)

        # Calculate H_k
        H_k = calculateBFGS_H(H_k, x_vec_current, x_vec_previous, g)
        # Calculate the descent direction
        p_vec = -H_k*g(x_vec_current)
        # Calculate the step length
        alpha = backtrackingLineSearch(object_function, x_vec_current, p_vec)

        x_vec_previous = x_vec_current
        x_vec_current = x_vec_current + alpha*p_vec
        diff_value = norm(g(x_vec))
        if diff_value < tolerance
            break
        end
    end
    return x_vec_current
end

function calculateBFGS_H(H_k_pre, x_vec::Vector, x_vec_pre::Vector, g::Function)
    s_k = x_vec - x_vec_pre
    y_k = g(x_vec) - g(x_vec_pre)
    I_matrix = Diagonal(ones(size(x_vec)[1]))
    if (Transpose(y_k)*s_k) == 0
        return H_k_pre
    end
    tho_k = 1 / (Transpose(y_k)*s_k)
    H_k = (I_matrix-tho_k*s_k*Transpose(y_k)) * H_k_pre * (I_matrix-tho_k*y_k*Transpose(s_k)) + tho_k*s_k*Transpose(s_k)
    return H_k
end