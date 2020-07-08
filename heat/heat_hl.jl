using LinearAlgebra
using Test
using SparseArrays

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

"""
furnish a triplet corresponding to the iterator F for the state equation
 a initial value `u0`, and the timestep
"""
function get_F(P::Int64)
    n = convert(Int64,(P-4) / 4)
    M = n * n
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

"""
computes the final value in time T for the solution
"""

function U_final(p::Array{Float64,1}, T::Float64)
    F, u, dt = get_F(size(p,1))
    @debug "dt = $dt, nb_it = ", floor(Int64,T/dt)
    for t=dt:dt:T
        u = F(u,p)
    end
    return  u
end

"""
returns the i-th vector of the canonical base of Rⁿ
"""
function 𝓔(n::Int64, i::Int64)
   return [1. * (k==i) for k=1:n]
end  

"""
computes the derivative of `F` function respect to `p` 
in the state equation 
Uⁿ⁺¹ = F(Uⁿ, p)
"""
function ∂ₚF_dense(p::Array{Float64,1})
    F, u0, dt = get_F(size(p, 1))
    P = size(p, 1)
    M = convert(Int64,((P-4)/4).^2)
    return [𝓔(M,i)'*F(zeros(M),𝓔(P,j)) for i=1:M, j=1:P]
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
    F, u0, dt = get_F(4n+4)
    return -dt .* (n+1).^2 * sparse(Is,Js,Vs, n*n,4n+4)
end 


@testset "check ∂ₚF" begin
    P = 16
    M = convert(Int64, ((P-4)/4.)^2)
    p = rand(P)
    F, u0, dt = get_F(P)
    U_rand = rand(M)
    ε = 1e-8
    fd_g = [𝓔(M,i)'*(F(U_rand,p+ ε * 𝓔(P,j))-F(U_rand,p)) / ε for i=1:M, j=1:P]
    @test norm(fd_g-∂ₚF(convert(Int64, √M))) < 10 * ε
    @test norm(fd_g-∂ₚF_dense(p)) < 10 * ε
end


"""
computes the partial differential of `F` function respect to `u`
"""
function ∂ᵤF_dense(M::Int64)
    P = 4 * convert(Int64, √M) + 4
    p = zeros(P)
    F, u0, dt = get_F(P)
    return [𝓔(M,i)'*F(𝓔(M,j),p) for i=1:M , j=1:M]
end 

"""
computes the partial differential of `F` function respect to `u` in a sparse way
"""
function ∂ᵤF(u::Array{Float64,1})
    M = size(u, 1)
    n = convert(Int64, √M)
    P = 4n + 4
    p = zeros(P)
    h = 1. /(n+1)
    dt = h.^2 /4
    U = ins⁻¹(u)
    set_boundary!(U,zeros(P))
    ΔU = heat_kernel(U) .* (n+1).^2 
    U[2:end-1,2:end-1] -= dt .* ΔU[2:end-1,2:end-1]
    return ins(U)
end 


@testset "check ∂ᵤF" begin
    P = 16
    M = convert(Int64, ((P-4)/4.)^2)
    p = rand(P)
    F, u0, dt = get_F(P)
    u = rand(M)
    fd_∇u = zeros(M,M)
    ε = 1e-8
    fd_∇u =[𝓔(M,i)'*(F(u+ε*𝓔(M,j),p)-F(u,p)) / ε for i=1:M ,j=1:M]
    @test norm(fd_∇u-∂ᵤF_dense(M)) < 10 .^2*ε
end

using Optim

function create_problem(target::Array{Float64,2}, T::Float64)
  n = size(target, 1)
  @assert size(target, 2) == n
  P = 4n+4
  M = n * n
  p0 = zeros(P)
  function simu(p::Array{Float64,1}, compute_gradient::Bool = true)
    u_f = U_final(p, T)
    v = sum((u_f - vec(target)).^2)
    λ = 2 .*(u_f - vec(target))
    if (!compute_gradient) 
        return v, λ
    end
    F, u, dt = get_F(P)
    fₚ = ∂ₚF(n)

    # academic aim 
    gₚ = zeros(P) 

    # just to learn the method , here u₀ does not depend on p 
    bₚ = zeros(M,P) 
    ∇f = zeros(P)

    for t=T:-dt:dt
        ∇f += fₚ'*λ
        # the adjoint equation is quasi the same a state equation
        λ = ∂ᵤF(λ) 
    end
    return v, gₚ + ∇f + bₚ' * λ
  end
  return simu, p0
end

@testset "check ∇g" begin

  n = 3
  target = rand(n, n)
  ε = 1e-8
  T = 10.
  simu, x0 = create_problem(target, T)
  P = size(x0, 1)
  f0, g0 = simu(x0, true) 
  g_fd = [(simu(x0+ε*𝓔(16,i), false)[1]-f0) / ε for i=1:P]
  @test norm(g_fd-g0)  < 100 * ε

end

using PyPlot
n = 15
T = 10.
x = [ 1. *i /(n+1) for i=1:n ]
y = copy(x)
XX = repeat(x,1,n)
YY = repeat(y',n,1)
Z = exp.(XX).*sin.(YY)
σ, p0 = create_problem(Z, T)
f = x-> σ(x, false)[1]

function ∇f!(g::Array{Float64,1}, p::Array{Float64,1})
    fp,gp = σ(p)
    g[:] = gp
end

res = optimize(f,∇f!, p0, LBFGS())
p_opt = Optim.minimizer(res)
u_opt = U_final(p_opt, T)
U = reshape(u_opt, n,n)

plot_surface(XX,YY,Z-U)
