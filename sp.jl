using NPBApps
using MPIClusterManagers
using Distributed

addprocs(MPIWorkerManager(4); threadleve=:multiple, mpiflags=`--map-by node --hostfile /home/ubuntu/hostfile`,  exeflags=`--threads=4 --check-bounds=no --optimize=3`)

@everywhere workers() using NPBApps
@everywhere workers() SP.go(SP.CLASS_W; timers=true)

