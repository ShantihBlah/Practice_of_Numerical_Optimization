include("element_func.jl")

abstract type AbstractNode end

mutable struct FunctionNode{T} <: AbstractNode
    value::T
    grad::T
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
        if i == length(forward_vec)
            node_vec[i].grad = 1.0
        else
            node_vec[i].grad = 0.0
        end
        for j in range(1, length=length(forward_vec[i]))
            push!(node_vec[i].forward, forward_vec[i][j])
        end
    end

    # println(vec_dim)

    return node_vec

end


function value_update(node_vec)
    for i in range(1, length=length(node_vec))
        if length(node_vec[i].backward) == 1
            node_vec[i].value = node_vec[i].f(node_vec[i].backward[1].value, node_vec[i].args)
        elseif length(node_vec[i].backward) == 2
            node_vec[i].value = node_vec[i].f(node_vec[i].backward[1].value, node_vec[i].backward[2].value, node_vec[i].args)
        end
    end
    return node_vec
end

function grad_update(node_vec)
    for i in range(1, length=length(node_vec))
        grad_index = length(node_vec) - i + 1
        if grad_index == length(node_vec)
            continue
        end
        expect_grad = 0.0
        for j in range(1, length=length(node_vec[grad_index].forward))
            next_node = node_vec[node_vec[grad_index].forward[j]]
            if length(next_node.backward) == 1
                expect_grad += next_node.grad * next_node.f_diff(node_vec[grad_index].value, next_node.args)
            else
                if node_vec[grad_index].value == next_node.backward[1].value
                    expect_grad += next_node.grad * next_node.f_diff(next_node.backward[1].value, next_node.backward[2].value, next_node.args)
                else
                    expect_grad += next_node.grad * next_node.f_diff(next_node.backward[2].value, next_node.backward[1].value, next_node.args)
                end
            end
        end
        node_vec[grad_index].grad = expect_grad
    end
    return node_vec
end


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

node_vec = value_update(node_vec)
node_vec = grad_update(node_vec)

println(node_vec[11].value)
grad_vec = [node_vec[1].grad; node_vec[2].grad]
println(grad_vec)
println(transpose(grad_vec)*p_vec)