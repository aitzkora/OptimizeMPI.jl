using MPI
MPI.Init()
println("Hi from $(MPI.Comm_rank(MPI.COMM_WORLD))!")
flush(stdout)
