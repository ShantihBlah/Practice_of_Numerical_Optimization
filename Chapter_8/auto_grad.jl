
abstract type AbstractNode end

abstract type LeafNode <: AbstractNode end

# mutable struct Variable{T} <: AbstractNode
#     value::T
#     grad::T

#     Variable(val::T) where T = new{T}(val, zero(val))
# end

mutable struct Variable{T} <: LeafNode
    value::T
    grad::T

    Variable(val::T) where T = new{T}(val)
    Variable(val::T, grad::T) where T = new{T}(val)
end

x = Variable([10.0; 4.0; 2.0])

println(x)

struct Node{FT <: Function, ArgsT <: Tuple, KwargsT <: NamedTuple} <: AbstractNode
    f::FT
    args::ArgsT
    kwargs::KwargsT
end

tupl3 = (a = 1, b = 2, c = "Hello Geeks")

node = Node(abs2, (1,2,3), tupl3)

println(node)

# # struct Node{FT <: Function, ArgsT <: Tuple, KwargsT <: NamedTuple} <: AbstractNode
# #     f::FT
# #     args::ArgsT
# #     kwargs::KwargsT
# # end

# abstract type Operator end

# module Trait
# import YAAD: Operator

# struct Method{FT} <: Operator
#     f::FT
# end

# struct Broadcasted{FT} <: Operator
#     f::FT
# end

# end # Trait

# struct Node{FT <: Operator, ArgsT <: Tuple, KwargsT <: NamedTuple} <: AbstractNode
#     f::FT
#     args::ArgsT
#     kwargs::KwargsT
# end

# # wrap function to Method
# Node(f::Function, args, kwargs) = Node(Trait.Method(f), args, kwargs)
# Node(op, args) = Node(op, args, NamedTuple())

# mutable struct CachedNode{NT <: AbstractNode, OutT} <: AbstractNode
#     node::NT
#     output::OutT
# end

# function CachedNode(f, args...; kwargs...)
#     node = Node(f, args, kwargs.data) # this constructs a Node
#     output = forward(node)
#     CachedNode(node, output)
# end

