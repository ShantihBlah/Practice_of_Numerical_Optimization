using AutoGrad

include("problems/unconstraint/quadratic_func.jl")
# include("problems/unconstraint/Rosenbrock_func.jl")

x = Param([10.0, 4.0, 2.0])

y = @diff sum(abs2, x)

println(value(y))
println(grad(y,x))