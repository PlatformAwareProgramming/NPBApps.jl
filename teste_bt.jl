using Distributed
process_count = [(2,16)#=,(3,4),(4,4),(5,4)=#]
addprocs(length(process_count))
@everywhere workers() using MPIClusterManagers
@everywhere workers() using MPI
for (w,np) in process_count
    fetch(@spawnat w addprocs(MPIWorkerManager(np); threadlevel=:multiple))
end
@everywhere workers() @everywhere workers() using NPBApps

using NPBApps
BT.go(BT.CLASS_D; itimer=2, npb_verbose=3#=, zone_mapping=[(2, [1,2,3]), (3, [4])]=#)

