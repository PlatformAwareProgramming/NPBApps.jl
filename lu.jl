using NPBApps
using MPIClusterManagers
using Distributed

addprocs(MPIWorkerManager(16); master_tcp_interface="172.31.1.0", threadleve=:multiple, mpiflags=`--map-by node --hostfile /home/ubuntu/hostfile`,  exeflags=`--threads=1 --check-bounds=no --optimize=3`)

@everywhere workers() using NPBApps
@everywhere workers() LU.go(LU.CLASS_C; timers=true)

