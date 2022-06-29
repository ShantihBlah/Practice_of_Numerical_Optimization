using LinearAlgebra
using ForwardDiff: gradient, hessian

include("backtracking_line_search.jl")


function lineSearch(object_function::Function, x_vec::Vector, max_iter_num::Int; BK_method="BFGS")
    """
    Line Search Method without contraints: 

    Paramerter ----
    object_function: the object function, i.e. f(x)
    x_vec: variant vector, i.e. x
    max_iter_num: the maximal number of iteration, i.e k-max
    BK_method: the methods of resolving B_k, ["Speedest_Descend", "Newtown", "BFGS"], "BFGS" is default value.

    Return ----
    x_vec_current: the optimal x, i.e. x*

    """
    tolerance = 1e-16
    x_vec_current = x_vec
    B_0 = Diagonal(ones(size(x_vec)[1]))
    B_k = B_0
    x_vec_previous = x_vec_current

    for i in range(1, length=max_iter_num)

        # Calculate B_k
        if (BK_method=="BFGS")
            B_k = BFGS(B_k, x_vec_current, x_vec_previous, object_function) # BFGS
        elseif (BK_method=="Newtown")
            B_k = hessian(object_function, x_vec_current) # the Newtown method
        elseif (BK_method=="Speedest_Descend")
            B_k = steepestDescentBkMatrix(x_vec_current) # the steepest descent method
        end

        # Calculate the descent direction
        p_vec = -inv(B_k)*gradient(object_function, x_vec_current)
        # Calculate the step length
        alpha = backtrackingLineSearch(object_function, x_vec_current, p_vec)
        
        x_vec_previous = x_vec_current
        x_vec_current = x_vec_current + alpha*p_vec
        diff_value = norm(gradient(object_function, x_vec))
        if diff_value < tolerance
            break
        end
    end
    return x_vec_current
end

function steepestDescentBkMatrix(x_vec::Vector)
    temp_vec = ones(size(x_vec)[1])
    return Diagonal(temp_vec)
end

function BFGS(B_k_pre, x_vec::Vector, x_vec_pre::Vector, object_function::Function)
    s_k = x_vec - x_vec_pre
    y_k = gradient(object_function, x_vec) - gradient(object_function, x_vec_pre)
    tolerance = 1e-16
    if norm(s_k) > tolerance
        B_k = B_k_pre - (B_k_pre*s_k*Transpose(s_k)*B_k_pre)/(Transpose(s_k)*B_k_pre*s_k) + (y_k*Transpose(y_k))/(Transpose(y_k)*s_k)
    else
        B_k = B_k_pre
    end
    return B_k
end
