using Combinatorics
using Plots
function p(m, n)
    if (m, n) == (1, 1)
        return 1/3
    end
    if n == 1
        return p(m-1, m-1)
    end
    if n == 0
        return 1
    end
    return sum(p(m, n-k) * (1/3)^k for k=1:n) + (2/3)^n * (1/2)
end
F(n) = 1 - p(n, n)
plot(f.(1:4))

function P(M, N)
    A = zeros(M, N)
    for m = 1:M, n=1:N
        if (m, n) == (1, 1)
            A[m, n] = 1/3
            continue
        end
        if n == 1
            A[m, n] = A[m-1, m-1]
            continue
        end
        A[m, n] = (1/2) * (2/3)^n + (1/3)A[m, n-1]*n
    end
    A[M, N]
end

function P(M, N)
    A = zeros(M, N)
    for m = 1:M, n=1:N
        if (m, n) == (1, 1)
            A[m, n] = 1/3
            continue
        end
        if n == 1
            A[m, n] = A[m-1, m-1]
            continue
        end
        A[m, n] = sum(A[m, n-k] * (1/3)^k for k=1:n-1) + (2/3)^n * (1/2) + (1/3)^n
    end
    A[M, N]
end

F(n) = 1 - P(n, n)
F(20)

plot([F.(1:10) f.(1:10)])

function P(m, n, N)
    λ = 1/2
    if N == 1
        return 1/3
    end
    if m + n == N
        n == 1 && return P(0, 0, N-1) * 1/3^m * 2/3^n
        return λ^m * (1 - λ)^n * 1/2
    end
    return P(m+1, n, N) + P(m, n+1, N)
end

F(n) = 1 - P(0, 0, n)
plot(F.(1:10))

P(0, 0, 2) == P(0, 1, 2) + P(1, 0, 2)
P(1, 0, 2) == P(2, 0, 2) + P(1, 1, 2)
