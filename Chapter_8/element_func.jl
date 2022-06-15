function power_func(x, args::Tuple)
    return x^args[1]
end

function power_diff_func(x, args::Tuple)
    return args[1]*(x^(args[1]-1))
end

function multiply_func(x, y, args::Tuple)
    return x*y   
end

function multiply_diff_func(x, y, args::Tuple)
    return y
end

function add_func(x, args::Tuple)
    return x+args[1]
end

function add_diff_func(x, args::Tuple)
    return 1
end

add_func(x, y, args::Tuple)=x+y

add_diff_func(x, y, args::Tuple) = 1
    
function sacle_func(x, args::Tuple)
    return args[1]*x
end

function sacle_diff_func(x, args::Tuple)
    return args[1]
end
