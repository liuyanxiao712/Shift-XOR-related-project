using Distributed
# addprocs(3)
@everywhere using AbstractAlgebra
@everywhere using Combinatorics
@everywhere using LinearAlgebra
@everywhere using BenchmarkTools
@everywhere using ProgressMeter

@everywhere T = 5
@everywhere shape = (12, 20)
@everywhere row, col = shape
@everywhere S, x = GF(2)["x"]
@everywhere φ = x^T - 1
@everywhere R = ResidueRing(S, φ)

@everywhere FF = GF(T)

@everywhere struct CircMat{T} <: MatElem{T}
    mat::MatElem{T}
    res::T
end

@everywhere function CircMat(mat, res=φ)
    new(mat, φ)
end

@everywhere function Generic.nrows(C::CircMat)
    return nrows(C.mat)
end

@everywhere function Generic.ncols(C::CircMat)
    return ncols(C.mat)
end

@everywhere function Base.getindex(C::CircMat, i::Int64, j::Int64)
    return getindex(C.mat, i, j)
end

@everywhere function Base.getindex(C::CircMat, ::Colon, a::Array{Int64, 1})
    return CircMat(getindex(C.mat, :, a), C.res)
end

@everywhere function Base.getindex(C::CircMat, ::Colon, a::Int)
    return CircMat(getindex(C.mat, :, a), C.res)
end

@everywhere function Base.similar(C::CircMat, i::Int64, j::Int64)
    return CircMat(similar(C.mat, i, j), C.res)
end

@everywhere function AbstractAlgebra.det(C::CircMat)
    return det(C.mat) % φ
end

@everywhere function issingular!(X::MatElem)
    if size(X, 1) ≠ size(X, 2)
        error("Not a square matrix")
    end

    row = size(X, 1)
    for i = 1:row

        # if not pivot
        if X[i, i] == 0
            #find next pivot
            found = false
            for k = i+1:row
                if X[k, i] ≠ 0
                    swap_rows!(X, i, k)
                    found = true
                end
            end
            found || return false
        end

        pivot = X[i, i]

        #row operations
        for j = i+1:row
            if X[j, i] ≠ 0
                exterior = X[j, i]
                for k = i:row
                    X[j, k] = X[i, k] * exterior - X[j, k] * pivot
                end
            end
        end
    end
    return true
end



@everywhere function randres(row=row, col=col, T=T)
    A = rand(1:T-1, row, col)
end



@everywhere function randcirc(row=row, col=col, T=T)
    A = rand(1:T-1, row, col)
    return CircMat(matrix(S, x.^A))
end

@everywhere function randmonoshift(row=row, col=col, T=T)
    A = rand(1:T-1, row, col)
    return matrix(S, x.^A)
end

@everywhere function randpolyshift(row=row, col=col, T=T)
    A = rand(0:2^(T) - 1, row, col)
    A = string.(A, base=2)
    strtopoly(str) = sum([str[i] == '0' ? 0 : x^(length(str)-i) for i = 1:length(str)])
    matrix(S, strtopoly.(A))
end

@everywhere function randmat(row=row, col=col, T=T)
    return matrix(GF(T), rand(1:T-1, (row, col)))
end


@everywhere function Base.getindex(X::MatElem{T}, c::Colon, a::Array{Int64,1}) where T <: AbstractAlgebra.RingElem
    S = [X[:, i:i] for i in a]
    A = S[1]
    for i = 2:size(S, 1)
        A = hcat(A, S[i])
    end
    return A
end

@everywhere function Base.getindex(X::MatElem{T}, a::Array{Int64, 1}, c::Colon) where T <: AbstractAlgebra.RingElem
    S = [X[i:i, :] for i in a]
    A = S[1]
    for i = 2:size(S, 1)
        A = vcat(A, S[i])
    end
    return A
end

@everywhere function Base.getindex(X::MatElem{T}, a::Array{Int64, 1}, b::Array{Int64, 1}) where T <: AbstractAlgebra.RingElem
    X[a, :][:, b]
end


@everywhere function Base.getindex(X::MatElem{T}, c::Colon, a::Int64) where T <: AbstractAlgebra.RingElem
    return getindex(X, :, [a]), X.res
end



@everywhere function isMDS(X)
    row, col = size(X)
    for comb in combinations(1:col, row)
        det(X[:, comb]) == 0 && return false
    end
    return true
end

@everywhere function isfullrank(X)
    row, col = size(X)
    for comb in combinations(1:col, row)
        det(X[:, comb]) ≠ 0 && return true
    end
    return  false
end

@everywhere function budget(X)
    row, col = size(X)
    for sub = row:col
        isfullrank(X[:, [1:sub...]]) && return sub
    end
    return col + 1
end

@everywhere function benchmark(cnt, total)
    println("$cnt/$total = $(cnt/total)")
end

@everywhere function MDStest(epoch, rand; shape=(row, col), T=T)
    MDS = 0
    MDS = @distributed (+) for i = 1:epoch
        isMDS(rand(row, col, T)) && 1
    end
    benchmark(MDS, epoch)
    MDS
end

@everywhere function budgettest(epoch, rand, shape=(row, col), T=T)
    row, col = shaperan
    votes = @distributed (+) for i = 1:epoch
        stat = Base.zeros(Int, col + 1 - row + 1)
        stat[budget(rand(shape..., T)) - row + 1] = 1
        stat
    end
    votes = Dict((((r + row - 1 == col + 1) ? 0 : r + row - 1) => (votes[r], votes[r]/epoch)) for r = 1:col-row + 2)
    votes = sort(Dict)
    votes
end

function main(args)
    rands = Dict("-c" => randcirc, "-s" => randpolyshift, "-ms" => randmonoshift)
    rand, epoch, row, col, T = args
    rand = rands[rand]
    epoch, row, col, T = parse.(Int, [epoch, row, col, T])

    @time budgettest(epoch, rand, (row, col), T)
end
# main(ARGS)

preal(n) = sum(isfullrank(randpolyshift(n, n, 2)) for i=1:1000) / 1000
#
elems = [0, 1, x, x+1]
#
reals = [preal.(1:20)...]


@everywhere function lexi(n)
    if n == 1
        return[[elem] for elem in elems]
    end


    combs = lexi(n-1)
    reshape([[elem comb] for elem in elems, comb in combs], 4^n)
end

# pfullrank(n) = sum(isfullrank(matrix(S, reshape(mat, n, n))) for mat in lexi(n^2)) / 4^(n^2)

function pfullrank(n)
    cnt = @distributed (+) for mat in lexi(n^2)
        isfullrank(matrix(S, reshape(mat, n, n)))
    end
    cnt / 4^(n^2)
end

# A = randpolyshift(4, 6, 2)
# isMDS(A)
# transpose(A)
#
# combs(n) = 4^(n^2)
# # t3 = @time pfullrank(3)
#
# expectedtime(n, t3) = combs(n) / combs(3)
# A = [1 1 1 1; x+1 0 x x+1; x+1 x+1 0 0; 0 0 x 1; 1 x+1 1 1; x+1 x x+1 1]


function adjugate(A)
    if size(A, 1) ≠ size(A, 2)
        error("Not a square matrix")
    end
    n = size(A, 1)
    C = matrix(S, zeros(Int, (n, n)))
    for i = 1:n, j=1:n
        C[i, j] = det(A[setdiff(1:n, i), setdiff(1:n, j)])
    end
    C
end


function randdiagshift(n)
    elems = [0, 1, x]
    A = matrix(S, zeros(Int, (n, n)))
    for i=1:n, j=1:n
        if i == j
            A[i, j] = x
        else
            A[i, j] = rand(0:1)
        end
    end
    A
end

det(randdiagshift(3 ))
