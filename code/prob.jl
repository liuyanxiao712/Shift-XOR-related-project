using Plots

reals = [ 0.746
 0.743
 0.798
 0.848
 0.9
 0.941
 0.965
 0.973
 0.983
 0.989]

# λ = 1.3
μ = 0
σ = 0.3

function p(k)
    k == 1 && return 0.25
    k == 2 && return 0.25
    k == 3 && return 0.2
    # ps = (1 - 3/4 * (1 - p(k-1))) * 1/4 + 3/4 * 1/2
    # ps = 1/4 + 3/4 * 1/2
    ps = σ
    ((1 + (2ps - 1)^(k-μ))/2)^(k-λ)
end

q(k) = 1 - p(k)

graph(n) = plot([q.(1:n) preal.(1:n)], label=["predicted", "real"])


function n(k, N)
    A = zeros(N)
    A[1] = 1
    for i = 2:N
        A[i] = 4^i
    end
    A[k]
end


pf(k, n) = prod(1 - 2.0^(-n + i) for i=0:k-1)


# EG(k, n) = pf(k, n) * 2^k + (1 - pf(k, n)) * 4^k
#
# P(k, n) = k == 1 ? (1 - 1/4^n) : prod(1 - EG(k, n)/4.0^n for k=1:n-1) * (1 - 1/4^n)
#
# plot([(n -> P(n, n)).(1:10), reals])
#
# plot((n -> pf(n, 2n)).(1:10))

# p(n) = 1 - (1/(2^n - 1))
#
# plot([p.(1:10), reals[1:end]])

using AbstractAlgebra
using Nemo

function randffmat(n)
    F, x = FiniteField(2, n+1, "x")
    S = [0, 1, x, x+1]
    S = F.(S)
    MS = MatrixSpace(F, n, n)
    M = MS(rand(S, (n, n)))
end

function statffmat(n, epochs)
    sum(det(randffmat(n)) == 0 ? 0 : 1 for i=1:epochs) / epochs
end

p(n) = statffmat(n, 1000)

plot([p.(1:10), reals])
