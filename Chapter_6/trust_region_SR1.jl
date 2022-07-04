using LinearAlgebra
using ForwardDiff: gradient, hessian

function trustRegionSR1(object_function::Function, x_vec::Vector, max_iter_num::Int; delta_max::Int = 10)
    tolerance = 1e-16
    delta_0 = 0.5 * delta_max
    eta = 0.5 * (10e-3)
    r = 0.5
    delta_k = delta_0
    zero_vec = zeros(size(x_vec)[1])

    # setup the 1st and 2nd Gradient
    f = object_function
    g = x_vec -> gradient(object_function, x_vec)
    H = x_vec -> hessian(object_function, x_vec)

    B_k = H(x_vec)

    for i in range(1, length=max_iter_num)

        # f_k = object_function(x_vec)
        # g_k = g(x_vec)
        # B_k = H(x_vec)

        # s_k = cauchyPointCalculation(g_k, delta_k, B_k)
        s_k = dogleg(g(x_vec), B_k, delta_k)

        y_k = g(x_vec+s_k) - g(x_vec)
        ared = f(x_vec) - f(x_vec+s_k)
        pred = -(Transpose(g(x_vec))*s_k + 0.5*Transpose(s_k)*B_k*s_k)
        
        if ared/pred > eta
            x_vec = x_vec + s_k
        else
            x_vec = x_vec
        end

        if ared/pred > 0.75
            if norm(s_k) <= 0.8*delta_k
                delta_k = delta_k
            else
                delta_k = 2.0 * delta_k
            end
        elseif ared/pred >= 0.1 && ared/pred <= 0.75
            delta_k = delta_k
        else
            delta_k = 0.5 * delta_k
        end
        if Transpose(s_k)*(y_k-B_k*s_k) >= r*norm(s_k)*norm(y_k-B_k*s_k)
            B_k = B_k + ((y_k-B_k*s_k)*Transpose(y_k-B_k*s_k)) / (Transpose(y_k-B_k*s_k)*s_k)
        else
            B_k = B_k
        end
    end
    return x_vec
end