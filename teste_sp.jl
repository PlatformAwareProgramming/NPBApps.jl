using Distributed
addprocs(2)
process_count = [(2,4),(3,4)#=,(4,4),(5,4)=#]
@everywhere workers() using MPIClusterManagers
@everywhere workers() using MPI
for (w,np) in process_count
    fetch(@spawnat w addprocs(MPIWorkerManager(np); threadlevel=:multiple))
end
@everywhere workers() @everywhere workers() using NPBApps

using NPBApps
SP.go(SP.CLASS_A; itimer=2, npb_verbose=3#=, zone_mapping=[(2, [1,3]), (3, [2,4])]=#)
