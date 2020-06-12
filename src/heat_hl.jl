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
    X = reshape(x,3,3)
    @test ins(ins⁻¹(x)) == x
    @test ins⁻¹(ins(X)) == [0 0 0; 0 X[2,2] 0; 0 0 0 ]
end


function set_boundary!(X::Array{Float64,2},g::Array{Float64,1})
    m = size(g,1)
    n = size(X,1)-2
    @assert size(X,1) == (m-4)/4 + 2
    @assert size(X,1) == size(X,2)
    X[1,:]         = g[1:n+2]
    X[2:end-1,1]   = g[n+3:2:3n+2] 
    X[2:end-1,end] = g[n+4:2:3n+2] 
    X[end,:]       = g[3n+3:4n+4]
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

@testset "boundary" begin
    a = [ 1 2 3; 4 0 5.; 6 7 8]
    @test get_boundary(a) == [1:8;]
    b = zeros(3,3)
    set_boundary!(b, [1:8.;])
    @test b ≈ a
end

function apply_Δ(x::Array{Float64,1}, g::Array{Float64,1})
    n = convert(Int64, √size(x,1))
    @assert size(g,1) == 4n + 4
    X = zeros(n+2, n+2)
    X[2:end-1,2:end-1] = reshape(x, n, n)
    set_boundary!(X, g)
    Y = heat_kernel(X) .* (n+1).^2
    return copy(ins(Y))
end

@testset "apply Δ" begin
    n = 4
    x = LinRange(0,1,n+2)
    xx = reshape(repeat(x, inner=n+2), n+2, n+2)
    yy = reshape(repeat(x, outer=n+2), n+2, n+2)
    u = yy.*(1 .- yy).*xx.^3
    g = get_boundary(u)
    Δu =  -6*xx.*yy.*(1 .-yy)+2*xx.^3
    @test norm(apply_Δ(ins(u),g) - ins(Δu) ) <= 1. /(n+1).^2
end

function compute(g::Array{Float64,1}, T::Float64)
    m = size(g,1)
    n = convert(Int64,(m-4) / 4)
    h = 1. /(n+1)
    dt = h.^2 /4
    u = zeros(n+2,n+2)
    set_boundary!(u, g)
    #println("dt = $dt, nb_it = ", floor(Int64,T/dt))
    for t=dt:dt:T
        Δu = heat_kernel(u) .* (n+1).^2
        u[2:end-1,2:end-1] -= dt .* Δu[2:end-1,2:end-1]
    end
    return u
end

using Optim
n = 5
u_cible = rand(n,n)
T = 10.
function simu(g::Array{Float64,1})
    u_calc = compute(g,T)
    return sum((u_cible[2:end-1,2:end-1] - u_calc[2:end-1,2:end-1]).^2)
end
g_0 = rand(4n-4)
res = optimize(simu, g_0)
sol = Optim.minimizer(res)
println("|u-u_cible|²=", simu(sol))
