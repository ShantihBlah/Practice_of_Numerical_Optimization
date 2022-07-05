using LinearAlgebra
using ForwardDiff: gradient, hessian

include("problems/unconstraint/quadratic_func.jl")
include("problems/unconstraint/rosenbrock_func.jl")
include("problems/unconstraint/linear_transform.jl")

include("Chapter_3/line_search.jl")
include("Chapter_4/trust_region.jl")
include("Chapter_5/linear_conjugate_gradient.jl")
include("Chapter_5/nonlinear_conjugate_gradient.jl")
include("Chapter_6/line_search_BFGS.jl")
include("Chapter_6/trust_region_SR1.jl")
include("Chapter_7/line_search_Newton_CG.jl")
include("Chapter_7/line_search_L_BFGS.jl")

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

function testNonlinearCG(problem, max_iter_num, algo_type, restart_flag)
    object_function = problem["obj_func"]
    init_vec = problem["init_vec"]
    init_y = object_function(init_vec)
    @time solution_vec = nonlinearConjugateGradient(object_function, init_vec, max_iter_num, algo_type=algo_type, restart_flag=restart_flag)
    solution_y = object_function(solution_vec)
    
    # println("Initial X: ", init_vec)
    println("Initial Y: ", init_y)
    # println("Final X: ", solution_vec)
    println("Final Y: ", solution_y)
end

function testBFGS(problem, max_iter_num)
    object_function = problem["obj_func"]
    init_vec = problem["init_vec"]
    init_y = object_function(init_vec)
    @time solution_vec = lineSearchBFGS(object_function, init_vec, max_iter_num)
    solution_y = object_function(solution_vec)
    
    # println("Initial X: ", init_vec)
    println("Initial Y: ", init_y)
    # println("Final X: ", solution_vec)
    println("Final Y: ", solution_y)
end

function testTrustRegionSR1(problem, max_iter_num)
    object_function = problem["obj_func"]
    init_vec = problem["init_vec"]
    init_y = object_function(init_vec)
    @time solution_vec = trustRegionSR1(object_function, init_vec, max_iter_num, delta_max=20)
    solution_y = object_function(solution_vec)
    
    # println("Initial X: ", init_vec)
    println("Initial Y: ", init_y)
    # println("Final X: ", solution_vec)
    println("Final Y: ", solution_y)
end

function testLineSearchNewtonCG(problem, max_iter_num)
    object_function = problem["obj_func"]
    init_vec = problem["init_vec"]
    init_y = object_function(init_vec)
    @time solution_vec = lineSearchNewtonCG(object_function, init_vec, max_iter_num, 1000)
    solution_y = object_function(solution_vec)
    
    # println("Initial X: ", init_vec)
    println("Initial Y: ", init_y)
    # println("Final X: ", solution_vec)
    println("Final Y: ", solution_y)
end

function testLBFGS(problem, max_iter_num)
    object_function = problem["obj_func"]
    init_vec = problem["init_vec"]
    init_y = object_function(init_vec)
    @time solution_vec = lineSearchLBFGS(object_function, init_vec, max_iter_num, 10)
    solution_y = object_function(solution_vec)
    
    # println("Initial X: ", init_vec)
    println("Initial Y: ", init_y)
    # println("Final X: ", solution_vec)
    println("Final Y: ", solution_y)
end

# Chapter_3
# testLineSearch(rosenbrock_setup, 100, "BFGS")

# Chapter_4
# testTrustRegion(rosenbrock_setup, 100)

# Chapter_5
# testPreliminaryCG(linear_transform_setup)
# testStandardCG(linear_transform_setup)
# testNonlinearCG(rosenbrock_setup, 1000, "Fletcher-Reeves", true)
# testNonlinearCG(rosenbrock_setup, 1000, "Hestenes-Stiefel", false)
# testNonlinearCG(rosenbrock_setup, 10000, "Polak-Ribiere", false)
# testNonlinearCG(rosenbrock_setup, 1000, "Dai-Yuan", true)

# Chapter_6
# testBFGS(rosenbrock_setup, 1000)
# testTrustRegionSR1(rosenbrock_setup, 10000)

# Chapter_7
# testLineSearchNewtonCG(rosenbrock_setup, 100)
testLBFGS(rosenbrock_setup, 1000)

