include("element_func.jl")

abstract type AbstractNode end

mutable struct FunctionNode{T} <: AbstractNode
    value::T
    Dp_value::T
    grad::Vector
    f::Function
    f_diff::Function
    args::Tuple
    forward::Vector
    backward::Tuple
    is_leaf::Bool
    FunctionNode(val::T) where T = new{T}(val, zero(val))
end


function initiate(func_vec, diff_func_vec, args_vec, backward_vec, forward_vec)
    node_vec = []
    vec_dim = 0
    for i in range(1, length=length(func_vec))
        x_node = FunctionNode(0.0)
        x_node.forward = []
        x_node.f = func_vec[i]
        x_node.args = args_vec[i]
        x_node.f_diff = diff_func_vec[i]
        if length(backward_vec[i]) == 1
            if backward_vec[i][1] == i
                x_node.backward = (x_node, )
                vec_dim = vec_dim + 1
                x_node.is_leaf = true
            else
                x_node.backward = (node_vec[backward_vec[i][1]],)
            end
        elseif length(backward_vec[i]) == 2
            x_node.backward = (node_vec[backward_vec[i][1]], node_vec[backward_vec[i][2]])
        end
        push!(node_vec, x_node)
    end

    for i in range(1, length=length(forward_vec))
        if node_vec[i].is_leaf == true
            node_vec[i].grad = zeros(vec_dim)
            node_vec[i].grad[i] = 1
        end
        for j in range(1, length=length(forward_vec[i]))
            push!(node_vec[i].forward, forward_vec[i][j])
        end
    end

    println(vec_dim)

    return node_vec

end


function update(node_vec)
    for i in range(1, length=length(node_vec))
        if length(node_vec[i].backward) == 1
            node_vec[i].value = node_vec[i].f(node_vec[i].backward[1].value, node_vec[i].args)
            node_vec[i].Dp_value = node_vec[i].f_diff(node_vec[i].backward[1].value, node_vec[i].args)*node_vec[i].backward[1].Dp_value
            node_vec[i].grad = node_vec[i].f_diff(node_vec[i].backward[1].value, node_vec[i].args)*node_vec[i].backward[1].grad
        elseif length(node_vec[i].backward) == 2
            node_vec[i].value = node_vec[i].f(node_vec[i].backward[1].value, node_vec[i].backward[2].value, node_vec[i].args)
            node_vec[i].Dp_value = node_vec[i].f_diff(node_vec[i].backward[1].value, node_vec[i].backward[2].value, node_vec[i].args)*node_vec[i].backward[1].Dp_value + 
                                   node_vec[i].f_diff(node_vec[i].backward[2].value, node_vec[i].backward[1].value, node_vec[i].args)*node_vec[i].backward[2].Dp_value
            node_vec[i].grad = node_vec[i].f_diff(node_vec[i].backward[1].value, node_vec[i].backward[2].value, node_vec[i].args)*node_vec[i].backward[1].grad + 
                               node_vec[i].f_diff(node_vec[i].backward[2].value, node_vec[i].backward[1].value, node_vec[i].args)*node_vec[i].backward[2].grad
        end
    end
    return node_vec
end


# # y = x1^2+x2^2+x3^2
# func_vec = [sacle_func, sacle_func, sacle_func, power_func, power_func, power_func, add_func, add_func]
# diff_func_vec = [sacle_diff_func, sacle_diff_func, sacle_diff_func, power_diff_func, power_diff_func, power_diff_func, add_diff_func, add_diff_func]
# args_vec = [(1.0,), (1.0,), (1.0,), (2.0,), (2.0,), (2.0,), (1.0,), (1.0,)]
# backward_vec = [(1,), (2,), (3,), (1,), (2,), (3,), (4, 5), (7, 6)]
# forward_vec = [(4,), (5,), (6,), (7,), (7,), (8,), (8,), (8,)]

# node_vec = initiate(func_vec, diff_func_vec, args_vec, backward_vec, forward_vec)


# node_vec[1].value = 1.
# node_vec[2].value = 1.
# node_vec[3].value = 1.


# p_vec = [3, 2, 2]
# node_vec[1].Dp_value = p_vec[1]
# node_vec[2].Dp_value = p_vec[2]
# node_vec[3].Dp_value = p_vec[3]

# node_vec = update(node_vec)
# println(node_vec[8].value, ", ", node_vec[8].Dp_value)
# println(node_vec[8].grad)


# Rosenbrock function
func_vec = [sacle_func, sacle_func, power_func, sacle_func, add_func, power_func, sacle_func, sacle_func, add_func, power_func, add_func]
diff_vec = [sacle_diff_func, sacle_diff_func, power_diff_func, sacle_diff_func, add_diff_func, power_diff_func, sacle_diff_func, sacle_diff_func, add_diff_func, power_diff_func, add_diff_func]
args_vec = [(1.0,), (1.0,), (2.0,), (-1.0,), (1.0,), (2.0,), (100.0,), (-1.0,), (1.0,), (2.0,), (1.0,)]
backward_vec = [(1,), (2,), (1,), (3,), (4,2), (5,), (6,), (1,), (8,), (9,), (10, 7)]
forward_vec = [(8,3), (5,), (4,), (5,), (6,), (7,), (11,), (9,), (10,), (11,), (11,)]

node_vec = initiate(func_vec, diff_vec, args_vec, backward_vec, forward_vec)

node_vec[1].value = 4.
node_vec[2].value = 2.


p_vec = [3; 2]
node_vec[1].Dp_value = p_vec[1]
node_vec[2].Dp_value = p_vec[2]

node_vec = update(node_vec)
println(node_vec[11].value, ", ", node_vec[11].Dp_value)
println(node_vec[11].grad)
println(transpose(node_vec[11].grad)*p_vec)