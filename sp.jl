using NPBApps
using MPIClusterManagers
using Distributed

addprocs(MPIWorkerManager(4); threadleve=:multiple)

@everywhere workers() using NPBApps
@everywhere workers() SP.go(SP.CLASS_W; timers=true)

