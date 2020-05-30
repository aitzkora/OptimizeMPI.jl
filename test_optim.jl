using Distributed
using Optim
using LinearAlgebra

function simu(n::Integer, x::Array{Float64,1})
    f = Ref{Float64}(0.)
    df = 0 
    c = cos.(1:n)[slice[1]:slice[2]]
    df = zeros(slice[2]-slice[1]+1)
    @assert  size(c,1) == size(x, 1)
    ccall((:compute_error, "./libpar_error.so"), 
          Cvoid, (Ref{Int32}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ptr{Float64}), size(x, 1), x, c, f, df)
    return f[], df
end

function partition(n::Integer, p::Integer)
     r = n % p
     m = ceil(Integer, n / p)
     part = fill( m, p )
     part[ (r+1) * (r> 0) + (p+1) * (r==0): p ] .= m - 1
     return part
end

function find_optim(n)

       comm = MPI.COMM_WORLD
       p = MPI.Comm_size(comm)
       r = MPI.Comm_rank(comm)
       global slice
       x = zeros(n)
       part = [0;cumsum(partition(n, p))] .+1
       slice = (part[r+1], part[r+2] -1)
       x[:] = 1:n
       x_loc = x[slice[1]:slice[2]];
       f,df = simu(n, x_loc)
       if (r == 0) 
           println("cost function at x0 = " , f)
       end
       cost = x-> simu(n,x)[1]
       grad = x-> simu(n,x)[2]
          
       res = optimize(cost, grad, x_loc, LBFGS(); inplace =false)
       x_min = Optim.minimizer(res)
       #x_glob = MPI.Gatherv(x_min, partition(n,p), 0, comm)
       #if (rank == 0)
       #    println(" |sol - cos(1:$n)| = ", norm(x_glob-cos.(1:n)))
       #end

end


