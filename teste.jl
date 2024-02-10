using Distributed
addprocs(2)
process_count = [(2,4),(3,4)]
@everywhere workers() using MPIClusterManagers
for (w,np) in process_count
    fetch(@spawnat w addprocs(MPIWorkerManager(np)))
end
@everywhere workers() @everywhere workers() using NPBApps

using NPBApps
SP.go(SP.CLASS_W; itimer=2, npb_verbose=3)

