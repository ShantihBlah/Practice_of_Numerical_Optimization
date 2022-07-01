using LinearAlgebra
using ForwardDiff: gradient, hessian

include("problems/unconstraint/quadratic_func.jl")
include("problems/unconstraint/Rosenbrock_func.jl")
include("problems/unconstraint/linear_transform.jl")

include("Chapter_3/line_search.jl")
include("Chapter_4/trust_region.jl")
include("Chapter_5/linear_conjugate_gradient.jl")
include("Chapter_5/nonlinear_conjugate_gradient.jl")

function testLineSearch(problem, max_iter_num, BK_method)
    object_function = problem["obj_func"]
    init_vec = problem["init_vec"]
    init_y = object_function(init_vec)
    @time solution_vec = lineSearch(object_function, init_vec, max_iter_num, BK_method=BK_method)
    solution_y = object_function(solution_vec)
    
    # println("Initial X: ", init_vec)
    println("Initial Y: ", init_y)
    # println("Final X: ", solution_vec)
    println("Final Y: ", solution_y)
end

function testTrustRegion(problem, max_iter_num)
    object_function = problem["obj_func"]
    init_vec = problem["init_vec"]
    init_y = object_function(init_vec)
    @time solution_vec = trustRegion(object_function, init_vec, max_iter_num, delta_max=20)
    solution_y = object_function(solution_vec)
    
    # println("Initial X: ", init_vec)
    println("Initial Y: ", init_y)
    # println("Final X: ", solution_vec)
    println("Final Y: ", solution_y)
end

function testPreliminaryCG(problem)
    A = problem["A"]
    x = problem["x"]
    b = problem["b"]
    @time solution_x = preliminaryConjugateGradient(A, x, b)
    println(0.5*transpose(solution_x)*A*solution_x)
end

function testStandardCG(problem)
    A = problem["A"]
    x = problem["x"]
    b = problem["b"]
    @time solution_x = standardConjugateGradient(A, x, b)
    println(0.5*transpose(solution_x)*A*solution_x)
end

function testNonlinearCG(problem, max_iter_num, algo_type)
    object_function = problem["obj_func"]
    init_vec = problem["init_vec"]
    init_y = object_function(init_vec)
    @time solution_vec = nonlinearConjugateGradient(object_function, init_vec, max_iter_num, algo_type=algo_type)
    solution_y = object_function(solution_vec)
    
    # println("Initial X: ", init_vec)
    println("Initial Y: ", init_y)
    # println("Final X: ", solution_vec)
    println("Final Y: ", solution_y)
end

# testLineSearch(rosenbrock_setup, 10000, "Speedest_Descend")
# testTrustRegion(quadratic_setup, 100)
# testPreliminaryCG(linear_transform_setup)
# testStandardCG(linear_transform_setup)
testNonlinearCG(rosenbrock_setup, 10000, "Polak-Ribiere")
