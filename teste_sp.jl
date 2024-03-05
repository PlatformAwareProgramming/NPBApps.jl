using Distributed
addprocs(2; lazy=false)
process_count = [(2,4),(3,4)#=,(4,4),(5,4)=#]
@everywhere workers() using MPIClusterManagers
@everywhere workers() using MPI
for (w,np) in process_count
    fetch(@spawnat w addprocs(MPIWorkerManager(np); threadlevel=:multiple, lazy=false))
end
@everywhere workers() @everywhere workers() using NPBApps

using NPBApps
SP.go(SP.CLASS_S; itimer=2, npb_verbose=3)
