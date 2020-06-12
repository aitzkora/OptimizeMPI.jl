using LinearAlgebra
using Test

function heat_kernel(X::Array{Float64,2})
    Y = copy(X)
    Y[2:end-1,2:end-1] .= 0.
    Y[2:end-1,2:end-1] = 4. .* X[2:end-1,2:end-1] -  
                               X[1:end-2,2:end-1] -
                               X[3:end  ,2:end-1] -
                               X[2:end-1,1:end-2] -
                               X[2:end-1,3:end]
    return Y
end

@testset "kernel heat" begin
    a = ones(5,5)
    a[2:end-1,2:end-1] .= 0.
    res = copy(a)
    res[2:end-1,2:end-1] = [ -2 -1 -2; -1 0 -1; -2 -1 -2 ]
    @test heat_kernel(a) ≈ res
end

function ins(X::Array{Float64,2})
    return X[2:end-1,2:end-1][:]
end


function ins⁻¹(x::Array{Float64,1})
    n = convert(Int64, √size(x,1))
    X=zeros(n+2,n+2)
    X[2:end-1,2:end-1]=x
    return X
end

@testset "ins"  begin
    x = rand(9)
    @test ins(ins⁻¹(x)) == X
end


function apply_Δ(x::Array{Float64,1}, g::Array{Float64,1})
    n = convert(Int64, √size(x,1))
    @assert size(g,1) == 4n + 4
    X = zeros(n+2, n+2)
    X[2:end-1,2:end-1] = reshape(x, n, n)
    X[1,:] = g[1:n+2]
    X[2:end-1,1] = g[n+3:2:3n+2]
    X[2:end-1,end] = g[n+4:2:3n+2]
    X[end,:] = g[3n+3:4n+4]
    Y = heat_kernel(X) / (n+1).^2
    return copy(ins(Y))
end

function get_boundary(X::Array{Float64,2})
    @assert size(X,1) > 2
    @assert size(X,1) == size(X,2)
    n = size(X,1) - 2
    g = zeros(4n+4)
    g[1:n+2] = X[1,:] 
    g[n+3:2:3n+2] = X[2:end-1,1] 
    g[n+4:2:3n+2] = X[2:end-1,end]
    g[3n+3:4n+4] = X[end,:] 
    return g
end

@testset "extract" begin
    a = [ 1 2 3; 4 0 5.; 6 7 8]
    @test get_boundary(a) == [1:8;]
end

@testset "apply Δ" begin
    n = 1000
    x = linrange(0,1,n+2)
    x = reshape(repeat(x, inner=n+2), n+2, n+2)
    y = copy(x)
    u = y.*(1 .- y).*x.^3
    g = get_boundary(u)
    Δu =  -6*x.*y.*(1 .-y)+2*x.^3
    @test apply_Δ(ins(u),g) ≈ ins(Δu)
end


