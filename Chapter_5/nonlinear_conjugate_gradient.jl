using LinearAlgebra
using ForwardDiff: gradient, hessian

function nonlinearConjugateGradient(object_function::Function, x_vec::Vector, max_iter_num::Int; algo_type="Fletcher-Reeves")
    tolerance = 1e-16

    # setup the 1st and 2nd Gradient
    g = x_vec -> gradient(object_function, x_vec)

    f_0 = object_function(x_vec)
    g_k = g(x_vec)
    p_k = -g_k
    x_vec_k = x_vec

    # Choose concrete algorithm
    if algo_type == "Fletcher-Reeves"
        calculateBeta = calculateFletcherReevesBeta
    elseif algo_type == "Polak-Ribiere"
        calculateBeta = calculatePolakRibiereBeta
    elseif algo_type == "Hestenes-Stiefel"
        calculateBeta = calculateHestenesStiefelBeta
    elseif algo_type == "Dai-Yuan"
        calculateBeta = calculateDaiYuanBeta
    elseif algo_type == "Hager-Zhang"
        calculateBeta = calculateHagerZhangBeta
    else
        calculateBeta = calculateFletcherReevesBeta
    end

    for i in range(1, length=max_iter_num)
        alpha_k = backtrackingLineSearch(object_function, x_vec_k, p_k)
        x_vec_k_1 = x_vec_k + alpha_k * p_k
        g_k_1 = g(x_vec_k_1)
        beta_k_1 = calculateBeta(g_k_1, g_k, p_k)
        p_k_1 = -g_k_1 + beta_k_1*p_k
        p_k = p_k_1
        x_vec_k = x_vec_k_1
        g_k = g_k_1
        diff_value = norm(g_k)
        if (diff_value < tolerance)
            break
        end
    end
    return x_vec_k
end

function calculateFletcherReevesBeta(g_k_1, g_k, p_k)
    return (Transpose(g_k_1) * g_k_1) / (Transpose(g_k) * g_k)
end

function calculatePolakRibiereBeta(g_k_1, g_k, p_k)
    return (Transpose(g_k_1) * (g_k_1-g_k)) / (Transpose(g_k_1) * g_k_1)
end

function calculateHestenesStiefelBeta(g_k_1, g_k, p_k)
    if (Transpose(g_k_1-g_k) * p_k) == 0
        return 0.0
    end
    return (Transpose(g_k_1) * (g_k_1-g_k)) / (Transpose(g_k_1-g_k) * p_k)
end

function calculateDaiYuanBeta(g_k_1, g_k, p_k)
    if (Transpose(g_k_1-g_k) * p_k) == 0
        return 0.0
    end
    return (Transpose(g_k_1) * g_k_1) / (Transpose(g_k_1-g_k) * p_k)
end

function calculateHagerZhangBeta(g_k_1, g_k, p_k)
    y_k = g_k_1 - g_k
    g_k_1_square = Transpose(g_k_1) * g_k_1
    denominator_term = Transpose(y_k) * p_k
    if denominator_term == 0
        return 0.0
    end
    return Transpose(y_k - 2*(g_k_1_square / denominator_term) * p_k) * (g_k_1 / denominator_term)
end

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