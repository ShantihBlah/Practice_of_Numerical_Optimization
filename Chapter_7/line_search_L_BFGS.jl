using LinearAlgebra
using ForwardDiff: gradient, hessian

include("../Chapter_3/backtracking_line_search.jl")

function lineSearchLBFGS(object_function::Function, x_vec::Vector, max_iter_num::Int, recent_iter_num::Int)
    tolerance = 1e-16
    x_vec_current = x_vec
    x_vec_previous = x_vec_current

    rollout = Vector{Vector}()
    gamma_k = 1.0

    # setup the 1st Gradient
    g = x_vec -> gradient(object_function, x_vec)

    for k in range(0, length=max_iter_num)
        H_0_k = gamma_k*Diagonal(ones(size(x_vec)[1]))
        # Calculate -H_k*g_k
        p_vec, gamma_k, rollout = calculateLBFGS_H(H_0_k, x_vec_current, x_vec_previous, g, recent_iter_num, rollout)
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

function calculateLBFGS_H(H_0_k, x_vec::Vector, x_vec_pre::Vector, g::Function, recent_iter_num::Int, rollout::Vector)
    s_k = x_vec - x_vec_pre
    y_k = g(x_vec) - g(x_vec_pre)
    if (Transpose(y_k)*s_k) == 0
        return -g(x_vec), 1, rollout
    end
    q = g(x_vec)
    m = length(rollout)
    for i in range(1, length=m)
        s_i = rollout[m+1-i][1]
        y_i = rollout[m+1-i][2]
        tho_i = 1 / (Transpose(y_i)*s_i)
        alpha_i = tho_i * Transpose(s_i) * q
        q = q - alpha_i*y_i
    end
    r = H_0_k*q
    for i in range(1, length=m)
        s_i = rollout[i][1]
        y_i = rollout[i][2]
        tho_i = 1 / (Transpose(y_i)*s_i)
        alpha_i = tho_i * Transpose(s_i) * q
        beta = tho_i*Transpose(y_i)*r
        r = r + (alpha_i-beta)*s_i
    end
    if length(rollout) < recent_iter_num
        push!(rollout, [s_k, y_k])
    else
        popfirst!(rollout)
        push!(rollout, [s_k, y_k])
    end
    gamma = (Transpose(s_k)*y_k) / (Transpose(y_k)*y_k)
    return -r, gamma, rollout
end