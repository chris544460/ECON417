using Distributions
using Interpolations
using Plots
include("gauher2.jl")

# Parameters
β = 0.7
r = 0.2
b = 1
σ = 0.5
n = 10

# Utility function
u(c) = log(c)

# Gaussian quadrature
weights, abscissas = gauher(n)

# Value function for t = 3
v3(s) = u(s)

# Decision rule for t = 2
function g2(s2)
    return (s2 * β * r) / (1 + β * r)
end

# Value function for t = 2
function v2(s2)
    a3 = g2(s2)
    return u(s2 - a3) + β * v3((1 + r) * a3)
end

# Bellman operator for t = 1
function bellman_operator(v, s1, a2)
    expectation = 0.0
    for i in 1:n
        u_i = σ * abscissas[i] * sqrt(2)
        y_next = b + exp(u_i)
        s2 = (1 + r) * a2 + y_next
        expectation += weights[i] * v(s2)
    end
    return sqrt(pi) * β * expectation
end

# Value function for t = 1 using cubic spline interpolation
function compute_v1(s1_values)
    v1_values = zeros(length(s1_values))
    for (i, s1) in enumerate(s1_values)
        objective(a2) = -u(s1 - a2) - bellman_operator(v2, s1, a2)
        v1_values[i] = -optimize(objective, 0, s1, GoldenSection()).minimum
    end
    return CubicSplineInterpolation(s1_values, v1_values)
end

# Set the grid for s1
s1_min = b
s1_max = b + exp(3 * σ) - 1
s1_values = range(s1_min, s1_max, length=100)

# Compute v1 using cubic spline interpolation
v1 = compute_v1(s1_values)

# Calculate v0(b) and g0(b)
function compute_v0_and_g0()
    objective(a1) = -u(b - a1) - bellman_operator(v1, b, a1)
    result = optimize(objective, 0, b, GoldenSection())
    v0_b = -result.minimum
    g0_b = result.minimizer
    return v0_b, g0_b
end

v0_b, g0_b = compute_v0_and_g0()

# Display the results
println("v0(b) = $v0_b")
println("g0(b) = $g0_b")

# Plot the value functions
plot(s1_values, v1.(s1_values), label="v1", xlabel="s1", ylabel="Value", legend=:topleft)
plot!(s1_values, v2.(s1_values), label="v2")
plot!(s1_values, v3.(s1_values), label="v3")

# Plot the decision rules
g1(s1) = optimize(a2 -> -u(s1 - a2) - bellman_operator(v1, s1, a2), 0, s1, GoldenSection()).minimizer
s_values = range(s1_min, s1_max, length=100)
plot(s_values, g1.(s_values), label="g1", xlabel="s", ylabel="Decision", legend=:topleft)
plot!(s_values, g2.(s_values), label="g2")

# Display the plots
display(plot)


   

