include("spline2.jl")
include("creategrid1.jl")
# include("gauss_quad.jl")
include("brent2.jl")
include("gauher2.jl")

using Plots, Random, Distributions
pyplot()

using LinearAlgebra

beta = 0.7
sigma = 0.5
mu = 0
r = 0.2
b = 1

function gaussian_quad(a_2, h, mu, sigma, n)
    w, x = gauher(n)
    sum = 0.0
 
    for i in 1:n
        s_2 = (1 + r) * a_2 + b + exp(sigma * x[i] * sqrt(2) + mu)
        sum += w[i] * h(s_2) / sqrt(pi)
    end
    return sum
end

function v_3(s_3)
    return log(s_3)
end

function v_2(s_2)
    return log(1 / (1 + beta) * s_2) + beta * log(beta / (1 + beta) * (1 + r) * s_2)
end

function g_2(s_2)
    return beta / (1 + beta) * s_2
end

function find_v1()
    N = 501
    S_1 = creategrid(b, (1 + r) * b + b + exp(3 * sigma), N)

    v_new = zeros(N)
    g = zeros(N)

    for i in 1:N
        v_new[i], g[i] = next_step(S_1[i])
    end

    return makespline(S_1, v_new), makespline(S_1, g)
end

function v_1(s_1)
    v1_spline = find_v1()[1]
    return interp(s_1, v1_spline)[1]
end

function v_0()
    N = 501
    s_0 = b

    v_0_max, a_1_max = brent(0, s_0 / 2, s_0, a_1 -> - (log(s_0 - a_1) + beta * gaussian_quad(a_1, s_1 -> v_1(s_1), mu, sigma, 10)), 1e-8, 1e-8)[1:2]
    return -v_0_max, a_1_max
end

function next_step(s_1)
    v_1_max, a_2_max = brent(0, s_1 / 2, s_1, a_2 -> - (log(s_1 - a_2) + beta * gaussian_quad(a_2, s_2 -> v_2(s_2), mu, sigma, 10)), 1e-8, 1e-8)[1:2]
    return -v_1_max, a_2_max
end

function compute_path()
    u = Normal(mu, sigma)

    t = [0, 1, 2, 3]
    s = zeros(4)
    v = zeros(4)

    s[1] = b
    v[1], a1 = v_0()

    v1_spline, g1_spline = find_v1()
    y1 = b + exp(rand(u, 1)[1])
    s[2] = (1 + r) * a1 + y1
    v[2] = interp(s[1], v1_spline)[1]
    a2 = interp(s[1], g1_spline)[1]

    s[3] = (1 + r) * a2 + b + exp(rand(u, 1)[1])
    v[3] = v_2(s[2])
    a3 = g_2(s[2])

    s[4] = (1 + r) * a3
    v[4] = v_3(s[3])

    # color red
    plot_v_path = plot(t, v, title = "Path of Value", label = "value", xlabel = "t", ylabel = "v", color = :red)
    plot_s_path = plot(t, s, title = "Path of Asset", label = "asset", xlabel = "t", ylabel = "s", color = :orange)

    display(plot(plot_s_path, plot_v_path, layout = (1, 2), size = (900, 600)))
    # save plot
    savefig("path.png")
end

function main()
    v1_spline, g1_spline = find_v1()
    v0, a1 = v_0()

    display(v0)
    display(a1)

    max_s1 = (1 + r) * b + b + exp(3 * sigma)
    nxfine = 200
    s1_vals = creategrid(b, max_s1, nxfine)
    g1 = broadcast(s1 -> interp(s1, g1_spline)[1], s1_vals)
    v1 = broadcast(s1 -> interp(s1, v1_spline)[1], s1_vals)

    max_s2 = (1 + r) * max_s1 + b + exp(3 * sigma)
    s2_vals = creategrid(b, max_s2, nxfine)
    v2 = broadcast(s2 -> v_2(s2), s2_vals)
    g2 = broadcast(s2 -> g_2(s2), s2_vals)

    max_s3 = (1 + r) * max_s2 + b + exp(3 * sigma)
    s3_vals = creategrid(b, max_s3, nxfine)
    v3 = broadcast(s3 -> v_3(s3), s3_vals)

    plot_g1 = plot(s1_vals, g1, title = "Optimal Savings Policy", label = "savings", xlabel = "s_1", ylabel = "savings", color = :red)
    plot_v1 = plot(s1_vals, v1, title = "Value Function", label = "v_1", xlabel = "s_1", ylabel = "value_func", color = :orange)
    plot_g2 = plot(s2_vals, g2, title = "Optimal Savings Policy", label = "a_3", xlabel = "s_2", ylabel = "a_3", color = :green)
    plot_v2 = plot(s2_vals, v2, title = "Value Function", label = "v_2", xlabel = "s_2", ylabel = "v_2", color = :blue)
    plot_v3 = plot(s3_vals, v3, title = "Value Function", label = "v_3", xlabel = "s_3", ylabel = "v_3", color = :purple)
    display(plot(plot_v1, plot_g1, plot_v2, plot_g2, plot_v3, layout = (3, 2), size = (800, 800)))
    # dwownload plot plot(plot_v1, plot_g1, plot_v2, plot_g2, plot_v3, layout = (3, 2), size = (800, 800)) as png

    savefig("plot.png")
end

# compute_path()

main()