using Distributed
using Optim
using LinearAlgebra

function create_par_problem(Z::Array{Float64,2}, T, proc)
  n = size(Z, 1)
  @assert size(Z, 2) == n
  P = 4n+4
  p0 = zeros(P)
  function cost(x::Array{Float64,1})
      f = Ref{Float64}(0.)
      @assert size(x) == (P,) 
      ccall((:compute_error, "./libheat_solve.so"), 
            Cvoid, (Ref{Int32}, Ref{Int32}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ptr{Float64}), 
            n,                proc,            T,            x,            Z,            f, convert(Ptr{Float64}, 0))
      return f[]
  end
  
  function grad!(g::Array{Float64,1}, x::Array{Float64,1})
      f = Ref{Float64}(0.)
      @assert size(x) == (P,) 
      @assert size(g) == (P,) 
      ccall((:compute_error, "./libheat_solve.so"), 
            Cvoid, (Ref{Int32}, Ref{Int32}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64},Ptr{Float64}), 
                    n,                proc,            T,            x,            Z,            f, g)
  end 
  return cost, grad! , p0
end

function find_optim(n)

       MPI.Init() 
       comm = MPI.COMM_WORLD
       p = MPI.Comm_size(comm)
       r = MPI.Comm_rank(comm)
    
       is_master = (r == 0)
       T = 1.
       x = [ 1. *i /(n+1) for i=1:n ]
       y = copy(x)
       XX = repeat(x,1,n)
       YY = repeat(y',n,1)
       Z = exp.(XX).*sin.(YY)

       σ, ∇σ!, x0 = create_par_problem(Z, T , convert(Int64, √p))
      
       # One evaluation of cost and gradient       
       f0 = σ(x0)
       g0 = zeros(size(x0))
       ∇σ!(g0, x0)
       if (is_master)
         println("f0 = $f0, |g0| = ", norm(g0))
       end
    
       #res = optimize(σ, ∇σ!, x0, ConjugateGradient(), Optim.Options(g_tol = 1e-7,
       #                            iterations = 1000, store_trace = false, show_trace = false ))
       #x_min = Optim.minimizer(res)
       #f_min = σ(x_min)
       #if (r == 0)
       #    println(" |f(x_min)| = ",  f_min)
       #end
end
