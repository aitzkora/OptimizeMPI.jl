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


"""
x = ins(X)

extract the interior of the matrix (X) and returns it as a vector
"""
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

function get_F(g::Array{Float64,1})
    m = size(g,1)
    n = convert(Int64,(m-4) / 4)
    h = 1. /(n+1)
    dt = h.^2 /4
    u0 = zeros(n*n)
    function F(u::Array{Float64,1},g::Array{Float64,1})
        U = ins⁻¹(u)
        set_boundary!(U,g)
        ΔU = heat_kernel(U) .* (n+1).^2 
        U[2:end-1,2:end-1] -= dt .* ΔU[2:end-1,2:end-1]
        return ins(U)
    end
    return F, u0, dt
end


function U_final(g::Array{Float64,1}, T::Float64)
    F, u, dt = get_F(g)
    @debug "dt = $dt, nb_it = ", floor(Int64,T/dt)
    for t=dt:dt:T
        u = F(u,g)
    end
    return  u
end

"""
computes the derivative of `F` function respect to `p` 
in the state equation 
Uⁿ⁺¹ = F(Uⁿ, p)
"""
function ∂ₚF_dense(p::Array{Float64,1})
    F,u0, dt = get_F(p)
    P = size(p, 1)
    M = convert(Int64,((P-4)/4).^2)
    g = zeros(M,P)
    for i=1:M
       for j=1:P
           g[i,j]= [1. * (k==i) for k=1:M]'*F(zeros(M),[1. * (k==j) for k=1:P]) 
       end
    end 
    return g
end 

"""
computes the constant (because F is linear) ∂ₚF in a sparse way
"""
function ∂ₚF(n::Int64)

    l2c = (i,j) -> 1 .+ (j .-2) .* n .+ (i .- 2)
    r1 = [3:n;]
    Is = [ l2c(2,r1) ; l2c(r1, n+1) ; l2c(r1,2.) ; l2c(n+1,r1) ]
    Js = [ r1 ; 2 .*r1 .+ n ; 2 .*r1 .+ (n-1) ; r1 .+ (3n+2)] 
    Vs =  -ones(4n-8)
    Is = [Is ; ones(2)  ; l2c(n+1,n+1)*ones(2); l2c(2,n+1)*ones(2) ; l2c(n+1,2)*ones(2) ] 
    Js = [Js ; [[2; n+3]; [3n+2 ; 4n+3 ]      ; [n+1; n+4]         ; [3n+1; 3n+4] ] ]
    Vs = [Vs ; -ones(8)]
    F, u0, dt = get_F(zeros(4*n+4))
    return -dt .* (n+1).^2 * sparse(Is,Js,Vs, n*n,4n+4)
end 


@testset "check ∂ₚF" begin
    P = 16
    M = convert(Int64, ((P-4)/4.)^2)
    p = rand(P)
    F, u0, dt = get_F(p)
    U_rand = rand(M)
    fd_g = zeros(M,P)
    ε = 1e-8
    for  i=1:M
       for  j=1:P
        fd_g[i,j]= [1. * (k==i) for k=1:M]'*(F(U_rand,p+ [ε * (k==j) for k=1:P])-F(U_rand,p)) / ε
        end
    end
    @test norm(fd_g-∂ₚF(convert(Int64, √M))) < 10 .^2*ε
    @test norm(fd_g-∂ₚF_dense(p)) < 10 .^2*ε
end

#using Optim
#n = 5
#u_cible = rand(n,n)
#T = 10.
function simu(g::Array{Float64,1})
    u_calc = U_final(g,T)
    return sum((u_cible[2:end-1,2:end-1] - u_calc[2:end-1,2:end-1]).^2)
end
#g_0 = rand(4n-4)
#res = optimize(simu, g_0)
#sol = Optim.minimizer(res)
#println("|u-u_cible|²=", simu(sol))
