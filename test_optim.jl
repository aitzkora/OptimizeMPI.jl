using Distributed
using Optim
using LinearAlgebra

function simu(x::Array{Float64,1})
    f = Ref{Float64}(0.)
    c = cos.(1:n)[slice[1]:slice[2]]
    df = zeros(slice[2]-slice[1]+1)
    @assert  size(c,1) == size(x, 1)
    ccall((:compute_error, "./libpar_error.so"), 
          Cvoid, (Ref{Int32}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ptr{Float64}), size(x, 1), x, c, f, df)
    return f[], df
end

function cost(x::Array{Float64,1})
    f,G = simu(x)
    return f
end

function grad(x::Array{Float64,1})
    f,G = simu(x)
    return G
end 

function find_optim(n)

       comm = MPI.COMM_WORLD
       nb_procs = MPI.Comm_size(comm)
       rank = MPI.Comm_rank(comm)
       global slice
       
       x = zeros(n)
       #slice = compute_slice(n, Int32(nb_procs), Int32(rank))
       println("je suis $rank", exp(rank))
       #x[:] = 1:n
       #x_loc = x[slice[1]:slice[2]]
       #f,df = simu(x_loc)
       #if (rank == 0) 
       #    println("cost function at x0 = " , f)
       #end
 
       #res = optimize(cost, grad, x_loc, LBFGS(); inplace =false)
       #x_min = Optim.minimizer(res)
       #x_glob = MPI.Gatherv(x_min, Cint[5,5], 0, comm)
       #if (rank == 0)
       #    println(" |sol - cos(1:$n)| = ", norm(x_glob-cos.(1:n)))
       #end

end


function compute_slice(n::Int32, procs::Int32, rank::Int32) 
    s = zeros(Int32, 2)
    length =  n / procs
    if ( mod(n, procs) == 0 ) 
        s[1] = rank * length + 1
        s[2] = s[1] + length - 1
    else
        if (rank < procs - 1)
            s[1] = rank * (length + 1) + 1
            s[2] = s[1] + (length + 1) - 1
        else
            s[1] = rank * (length + 1) + 1
            s[2] = s[1] + mod(n, length + 1) - 1
        end 
    end
    return s
end 
