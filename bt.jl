using NPBApps
using MPIClusterManagers
using Distributed

addprocs(MPIWorkerManager(1); #=master_tcp_interface="172.31.1.0", threadlevel=:multiple, mpiflags=`--map-by node --hostfile /home/ubuntu/hostfile`,  exeflags=`--threads=4 --check-bounds=no --optimize=3`=#)

@everywhere workers() using NPBApps
@everywhere workers() BT.go(BT.CLASS_B; timers=true)

