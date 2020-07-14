using Distributed
using Optim
using LinearAlgebra

function simu!(n::Integer, x::Array{Float64,1}, df::Array{Float64,1})
    f = Ref{Float64}(0.)
    c = cos.(1:n)[slice[1]:slice[2]]
    #@assert  size(c,1) == size(x, 1)
    ccall((:compute_error, "./libpar_error.so"), 
          Cvoid, (Ref{Int32}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ptr{Float64}), size(x, 1), x, c, f, df)
    return f[]
end

function partition(n::Integer, p::Integer)
  r = n % p
  m = ceil(Integer, n / p)
  part = fill( m, p )
  part[ (r+1) * (r> 0) + (p+1) * (r==0): p ] .= m - 1
  return part
end

function find_optim(n, script = true)

  # MPI preamble
  if (script)
    MPI.Init()
  end
  
  comm = MPI.COMM_WORLD
  p = MPI.Comm_size(comm)
  r = MPI.Comm_rank(comm)

  # computing partition and filling arrays
  global slice
  x = zeros(n)
  part = [0;cumsum(partition(n, p))] .+1
  slice = (part[r+1], part[r+2] -1)
  x[:] = zeros(n)
  x_loc = x[slice[1]:slice[2]];
  df = zeros(slice[2]-slice[1]+1)
  df_dummy = copy(df)
  
  # defini the simu! marginal functions
  cost = x -> simu!(n, x, df_dummy)
  grad! = (g,x)-> simu!(n,x,g)
  
  # evaluation of cost and gradient in x_0
  f = cost(x_loc)
  grad!(df, x_loc)

  # beware Cint here is mandatory
  df_glob = MPI.Gatherv(df, Cint[i for i in partition(n,p)], 0, comm)
  if (r == 0)
      println("f(x₀) = $f, |∇f(x₀)| = ", norm(df_glob))
  end 

  # optimization phas
  res = optimize(cost, grad!, x_loc, LBFGS(), Optim.Options(g_tol = 1e-12,
                               iterations = 1000, store_trace = false, show_trace = false ))

  # post optimization processing
  x_min = Optim.minimizer(res)
  x_glob = MPI.Gatherv(x_min,Cint[i for i in partition(n,p)] , 0, comm)

  if (r == 0)
      println("|sol - cos(1:$n)| = ", norm(x_glob-cos.(1:n)))
  end

  return res, x_glob
end
