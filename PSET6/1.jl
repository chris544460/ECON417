include("spline2.jl")


using Plots
using Statistics
# function to create the grid
function creategrid(a,b,n)
    return range(a, stop=b, length=n)
end

# function to compute the approximation error
function approx_error(f, fapprox)
    return abs.(f .- fapprox)
end

# function to perform linear interpolation
function linear_interp(x, y, xval)
    # perform linear interpolation without searchsortedfirst or searchsortedlast

    k = 1
    while xval > x[k+1]
        k += 1
    end
    return y[k] + (y[k+1] - y[k])/(x[k+1] - x[k])*(xval - x[k])
end

# function to compute the cubic spline
function cubic_spline(x, y)
    yspline = makespline(x, y)
    return x -> interp(x, yspline)[1]
end

# set up the grid and calculate the function values on the grid
N = 5
a = 0
b = 2π
x = creategrid(a, b, N)
f = sin.(x)
# try f = 2x +1
# f = 2x .+ 1
# perform linear interpolation
xval = range(a, stop=b, length=100)
fapprox_linear = broadcast(xval -> linear_interp(x, f, xval), xval) #linear_interp.(x, f, xval)

# broadcast(mm -> linear_interp(x, f, mm), xval)

# compute cubic spline interpolation
fapprox_cubic = cubic_spline(x, f).(xval)

# compute the approximation errors
error_linear = approx_error(sin.(xval), fapprox_linear)
error_cubic = approx_error(sin.(xval), fapprox_cubic)

# compute the average and maximum errors
avg_error_linear = mean(error_linear)
max_error_linear = maximum(error_linear)
avg_error_cubic = mean(error_cubic)
max_error_cubic = maximum(error_cubic)

# plot the results
plot(xval, error_linear, label="Linear Interpolation")
plot!(xval, error_cubic, label="Cubic Spline Interpolation")
xlabel!("x")
ylabel!("Approximation Error")
title!("Approximation Errors for sin(x)")

# now print them to the console
println("Average Error (Linear): $(round(avg_error_linear, digits=6))")
println("Max Error (Linear): $(round(max_error_linear, digits=6))")
println("Average Error (Cubic): $(round(avg_error_cubic, digits=6))")
println("Max Error (Cubic): $(round(max_error_cubic, digits=6))")

# plot the two functions
plot(xval, sin.(xval), label="sin(x)")
plot!(xval, fapprox_linear, label="Linear Interpolation")
plot!(xval, fapprox_cubic, label="Cubic Spline Interpolation")