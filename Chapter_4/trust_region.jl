using LinearAlgebra
using ForwardDiff: gradient, hessian

function trustRegion(object_function::Function, x_vec::Vector, max_iter_num::Int; delta_max = 10)
    tolerance = 1e-16
    delta_0 = 0.5 * delta_max
    eta = 0.5 * (1/4)
    delta_k = delta_0
    zero_vec = zeros(size(x_vec)[1])
    println(zero_vec)
    for i in range(1, length=max_iter_num)

        f_k = object_function(x_vec)
        g_k = gradient(object_function, x_vec)
        B_k = hessian(object_function, x_vec)

        # p_k = cauchyPointCalculation(g_k, delta_k, B_k)
        p_k = dogleg(g_k, B_k, delta_k)

        m_0 = resolveMk(f_k, g_k, B_k, zero_vec)
        m_k = resolveMk(f_k, g_k, B_k, p_k)
        rho_k = (object_function(x_vec)-object_function(x_vec+p_k)) / (m_0 - m_k)
        # println("rho_k: ", rho_k)

        if (rho_k < 0.25)
            delta_k = 0.25 * delta_k
        else
            if (rho_k > 0.75 && norm(p_k)==delta_k)
                delta_k = min(2*delta_k, delta_max)
            end
        end
        if (rho_k > eta)
            x_vec = x_vec+p_k
        end
        diff_value = norm(gradient(object_function, x_vec))
        if (diff_value < tolerance)
            break
        end
    end
    return x_vec
end

function cauchyPointCalculation(g_k::Vector, delta_k::Float64, B_k::Matrix)
    tau_k = 1
    sec_item = Transpose(g_k) * B_k * g_k
    if (sec_item > 0)
        tau_k = min(norm(g_k)^3 / (delta_k*sec_item), 1)
    end
    p_k = - tau_k * (delta_k/norm(g_k)) * g_k
    return p_k
end

function dogleg(g_k::Vector, B_k::Matrix, delta_k::Float64)
    sec_item = Transpose(g_k)*B_k*g_k
    p_U = - (Transpose(g_k)*g_k / sec_item) * g_k
    p_B = - inv(B_k) * g_k
    # Solve tau
    tau_k = 2.
    if (norm(p_B) > delta_k)
        # tau_k = 1.1
        # quad_no_neg: 0.25*(b^2-4*ac)
        quad_no_neg = (Transpose(p_B-p_U)*p_U)^2 - norm(p_B-p_U)^2*(norm(p_U)^2-delta_k^2)
        # tau should be in (1,2]
        if (quad_no_neg > 0)
            tau_k = (-Transpose(p_B-p_U)*p_U + sqrt(quad_no_neg)) / norm(p_B-p_U)^2 + 1
            if tau_k>2
                tau_k = 2
            end
        # tau should be in (0,1]
        else
            tau_k = delta_k/norm(p_U)
            if tau_k>1
                tau_k = 1
            end
        end
    end
    # println("tau_k: ", tau_k)
    # Solve p_tau
    p_tau = tau_k*p_U
    if (tau_k>1 && tau_k<=2)
        p_tau = p_U + (tau_k-1)*(p_B-p_U)
    end
    return p_tau
end

function resolveMk(f_k::Float64, g_k::Vector, B_k::Matrix, p_k::Vector)
    m_k = f_k + Transpose(g_k)*p_k + 0.5*Transpose(p_k)*B_k*p_k
    return m_k 
end

