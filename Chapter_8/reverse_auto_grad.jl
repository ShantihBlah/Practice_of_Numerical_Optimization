
abstract type LeafNode <: AbstractNode end


mutable struct Variable{T} <: LeafNode
    value::T
    grad::T

    Variable(val::T) where T = new{T}(val, zero(val))
    Variable(val::T, grad::T) where T = new{T}(val, grad)
end


struct Node{FT <: Function, DF <: Function, ArgsT <: Tuple} <: AbstractNode
    f::FT
    f_diff::DF
    args::ArgsT
end

function sacle_func(x, args::Tuple)
    return args[1]*x
end

function sacle_diff_func(x, args::Tuple)
    return args[1]
end

function test_func(x::Variable)
    node_1 = Node(sacle_func, sacle_diff_func, (3,))
    node_2 = Node(sacle_func, sacle_diff_func, (2,))
    node_3 = Node(sacle_func, sacle_diff_func, (2,))
    x.value[3] = node_1.f(x.value[1], node_1.args)+node_2.f(x.value[2], node_2.args)
    x.grad[1] = node_1.f_diff(x.value[1], node_1.args)
    x.grad[2] = node_2.f_diff(x.value[2], node_2.args)
    return x.value[3]
end

x = Variable([1.,2.,0.0], [1.,1.,1.])

y = test_func(x)

println(x)