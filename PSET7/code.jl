include("spline2.jl")
include("creategrid1.jl")
include("brent2.jl")
include("gauher2.jl")

# Collaborators: Amy, Kelvin

using Plots, Random, Distributions, LinearAlgebra

const BETA = 0.7
const SIGMA = 0.5
const MU = 0
const R = 0.2
const B = 1

function gaussian_quad(a_2, h, mu, sigma, n)
    w, x = gauher(n)
    return sum(w[i] * h((1 + R) * a_2 + B + exp(sigma * x[i] * sqrt(2) + mu)) / sqrt(pi) for i in 1:n)
end

function v_3(s_3)
    return log(s_3)
end

function v_2(s_2)
    return log(1 / (1 + BETA) * s_2) + BETA * log(BETA / (1 + BETA) * (1 + R) * s_2)
end

function g_2(s_2)
    return BETA / (1 + BETA) * s_2
end

function find_v1()
    N = 501
    S_1 = creategrid(B, (1 + R) * B + B + exp(3 * SIGMA), N)
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
    s_0 = B
    v_0_max, a_1_max = brent(0, s_0 / 2, s_0, a_1 -> -(log(s_0 - a_1) + BETA * gaussian_quad(a_1, s_1 -> v_1(s_1), MU, SIGMA, 10)), 1e-8, 1e-8)[1:2]
    return -v_0_max, a_1_max
end

function next_step(s_1)
    v_1_max, a_2_max = brent(0, s_1 / 2, s_1, a_2 -> -(log(s_1 - a_2) + BETA * gaussian_quad(a_2, s_2 -> v_2(s_2), MU, SIGMA, 10)), 1e-8, 1e-8)[1:2]
    return -v_1_max, a_2_max
end

function compute_path()
    u = Normal(MU, SIGMA)
    t = [0, 1, 2, 3]
    s = zeros(4)
    v = zeros(4)

    s[1] = B
    v[1], a1 = v_0()

    v1_spline, g1_spline = find_v1()
    y1 = B + exp(rand(u, 1)[1])
    s[2] = (1 + R) * a1 + y1
    v[2] = interp(s[1], v1_spline)[1]
    a2 = interp(s[1], g1_spline)[1]

    s[3] = (1 + R) * a2 + B + exp(rand(u, 1)[1])
    v[3] = v_2(s[2])
    a3 = g_2(s[2])

    s[4] = (1 + R) * a3
    v[4] = v_3(s[3])

    plot_v_path = plot(t, v, title="Path of Value", label="value", xlabel="t", ylabel="v", color=:red)
    plot_s_path = plot(t, s, title="Path of Asset", label="asset", xlabel="t", ylabel="s", color=:orange)

    display(plot(plot_s_path, plot_v_path, layout=(1, 2), size=(900, 600)))
    savefig("path.png")
end

function main()
    v1_spline, g1_spline = find_v1()
    v0, a1 = v_0()

    display(v0)
    display(a1)

    max_s1 = (1 + R) * B + B + exp(3 * SIGMA)
    nxfine = 200
    s1_vals = creategrid(B, max_s1, nxfine)
    g1 = broadcast(s1 -> interp(s1, g1_spline)[1], s1_vals)
    v1 = broadcast(s1 -> interp(s1, v1_spline)[1], s1_vals)

    max_s2 = (1 + R) * max_s1 + B + exp(3 * SIGMA)
    s2_vals = creategrid(B, max_s2, nxfine)
    v2 = broadcast(s2 -> v_2(s2), s2_vals)
    g2 = broadcast(s2 -> g_2(s2), s2_vals)

    max_s3 = (1 + R) * max_s2 + B + exp(3 * SIGMA)
    s3_vals = creategrid(B, max_s3, nxfine)
    v3 = broadcast(s3 -> v_3(s3), s3_vals)

    plot_g1 = plot(s1_vals, g1, title="Optimal Savings Policy", label="a2", xlabel="s1", ylabel="a2")
    plot_v1 = plot(s1_vals, v1, title="Value Function", label="v1", xlabel="s1", ylabel="v1")
    plot_g2 = plot(s2_vals, g2, title="Optimal Savings Policy", label="a3", xlabel="s2", ylabel="a3")
    plot_v2 = plot(s2_vals, v2, title="Value Function", label="v2", xlabel="s2", ylabel="v2")
    plot_v3 = plot(s3_vals, v3, title="Value Function", label="v3", xlabel="s3", ylabel="v3")
    display(plot(plot_v1, plot_g1, plot_v2, plot_g2, plot_v3, layout=(3, 2), size=(800, 800)))

    savefig("plot.png")
end

main()

compute_path()

