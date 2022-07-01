using ForwardDiff: GradientConfig, Chunk, gradient!, hessian, gradient, hessian
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

# x = rand(10000);
# out = similar(x);
# cfg1 = GradientConfig(rosenbrock, x, Chunk{1}());
# cfg4 = GradientConfig(rosenbrock, x, Chunk{4}());
# cfg10 = GradientConfig(rosenbrock, x, Chunk{10}());

# @time gradient!(out, rosenbrock, x, cfg1);
# @time gradient!(out, rosenbrock, x, cfg4);
# @time gradient!(out, rosenbrock, x, cfg10);


# x = [4.;2.];
# out = similar(x);
# cfg1 = GradientConfig(rosenbrock, x, Chunk{1}());
# @time gradient!(out, rosenbrock, x, cfg1);
# println(gradient!(out, rosenbrock, x, cfg1))
# hessian_matrix = hessian(rosenbrock, x)
# println(hessian_matrix)

function quadratic(vec)
    A = [1. 0. 2.; 1. 2. 0.; 5. 0. 2.]
    Q = A * Transpose(A)
    result = 0.5*transpose(vec)*Q*vec
    return result
end

# x = [10.;4.;2.];
x = [-9.0, 4.2, 8.]
println(quadratic(x))
out = similar(x);
cfg1 = GradientConfig(quadratic, x, Chunk{1}());
@time gradient!(out, quadratic, x, cfg1);
println(gradient!(out, quadratic, x, cfg1))
hessian_matrix = hessian(quadratic, x)
println(hessian_matrix)

g = x -> gradient(quadratic,x)
hessian_x = x -> hessian(quadratic, x)

x = [-9.0, 4.2, 8.]

println(g(x))
println(hessian_x(x))