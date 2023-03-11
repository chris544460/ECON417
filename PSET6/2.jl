using LinearAlgebra

include("spline2.jl")

# import optimal_value_function from 1a.jl in the same directory
function optimal_value_function(alpha, beta, delta, N, tolerance)
    k_ss = (alpha/(1/beta-1+delta))^(1/(1-alpha)) # steady-state capital stock
    k = range(0.5*k_ss, 1.5*k_ss, length=N) # grid for capital
    v = zeros(N) # initial set of values for the value function (step i)
    distance = Inf

    # optimal choice return value
    g_return = fill(0, N)

    while distance > tolerance
        v_new = zeros(N) # new set of values for the value function
        g = fill(0, N) # optimal choice

        for i in 1:N
            u = -Inf # initialize the value of the utility function (step ii)
            ki = k[i] # ki is the level of capital at the i-th grid point

            for j in 1:N
                kj = k[j] # kj is the level of capital at the j-th grid point
                if ki^alpha - kj + (1-delta)*ki>=0 # impose the constraint that consumption cannot be negative (step ii)
                    if log(ki^alpha - kj + (1-delta)*ki) + beta * v[j] > u # check if the utility function is increasing (step ii)
                        u = log(ki^alpha - kj + (1-delta)*ki) + beta * v[j] # calculate the value of the utility function (step ii)
                        g[i] = j # calculate the optimal choice (step iii)
                    end
                end
            end

            v_new[i] = u #+ beta * v[g[i]] # calculate a new set of values for the value function (step iii)
        end

        distance = maximum(abs.(v_new - v)) # compute the distance between the updated value function and the previous value function (step iv)
        v = v_new # update the value function (step iv)
        g_return = g # update the optimal choice
    end

    return v, g_return
end

using Interpolations

# Define the optimal decision rule using cubic spline interpolation
function q(k, g, k_grid)
    # make cubic spline interpolation
    g_spline = makespline(k_grid, g)
    # return the interpolated value
    return interp(k, g_spline)[1]
end

# Set parameter values
alpha = 0.4
beta = 0.9
delta = 0.1
tolerance = 1e-6

# Compute the steady-state capital stock
k_star = ((1/beta - 1 + delta)/alpha)^(1/(alpha-1)) # steady-state capital stock

# Set up the grid for capital
N = 5000
k_min = 1 # lower bound for the grid
k_max = 4 # upper bound for the grid
k_grid = range(k_min, stop=k_max, length=N) # grid for capital

# Find the index of the capital stock closest to k_star
index_k_star = searchsortedfirst(k_grid, k_star) # find the index of the capital stock closest to k_star

# Choose an initial value for k
k_init = 4
# Compute the optimal decision rule
v, g = optimal_value_function(alpha, beta, delta, N, tolerance) # compute the optimal value function and the optimal decision rule

# Compute the optimal dynamic path for the capital stock using cubic spline interpolation
k_min_new = 0.5 * k_star # lower bound for the new grid
k_max_new = 1.5 * k_star # upper bound for the new grid
N_new = 5000 # number of points on the new grid
k_grid_new = range(k_min_new, stop=k_max_new, length=N_new) # new grid for capital
k_path = zeros(N_new) # initialize the optimal dynamic path for the capital stock
k_path[1] = k_init # set the initial value of k

for i in 2:100
    k_path[i] = q(k_path[i-1], k_grid_new[g], k_grid_new) # compute the optimal dynamic path for the capital stock using the fitted cubic spline
end

# import plot from Plots.jl
using Plots

# Plot the optimal dynamic path for the capital stock (from k_path at time 1 to k_path at time 50
plot(k_path[1:100], title="Optimal dynamic path for the capital stock", xlabel="Time", ylabel="Capital stock", label="Capital stock", legend=:topleft)

# Confirm that the optimal dynamic path converges to the steady state
println("Distance between initial capital stock and steady state: ", abs(k_init - k_star))
println("Distance between final capital stock and steady state: ", abs(k_path[end] - k_star))
