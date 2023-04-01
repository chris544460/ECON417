include("gauher2.jl")

function expected_exponential(mu, sigma, n)
    (weights, abscissas) = gauher(n)
    g(z) = exp(z)
    result = 0.0

    for i in 1:n
        result += weights[i] * g(sigma * sqrt(2) * abscissas[i] + mu)
    end
    
    return sqrt(pi) * result
end

mu = 1
sigmas = [1, 2, 3]
n_values = 2:10

for sigma in sigmas
    println("Sigma: $sigma")
    for n in n_values
        expected_value = expected_exponential(mu, sigma, n)
        println("n: $n, Expected Value: $expected_value")
    end
    println()
end
