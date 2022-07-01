#
#  Rosenbrock function: https://en.wikipedia.org/wiki/Rosenbrock_function
#

using LinearAlgebra

function rosenbrock(x)
    a = one(eltype(x))
    b = 100 * a
    result = zero(eltype(x))
    for i in 1:length(x)-1
        result += (a - x[i])^2 + b*(x[i+1] - x[i]^2)^2
    end
    return result
end

# rosenbrock_init_vec = [4.;2.]
rosenbrock_init_vec = [4.;2.;2.;5.;8.;1.;3.;-1.;2.;-4.]
# rosenbrock_init_vec = rand(20)

rosenbrock_setup = Dict("obj_func" => rosenbrock, "init_vec" => rosenbrock_init_vec)