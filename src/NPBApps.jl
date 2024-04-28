module NPBApps

module SP

    using FortranFiles
    using OffsetArrays
    using Parameters
    using Printf
    using MPI   
    using StaticArrays
    using LoopVectorization
    using MPIClusterManagers
    using Distributed
    #using Traceur
    #using JET       
    #using InteractiveUtils
    #using ProfileView
    #using Semaphores

    include("common/DFVariable.jl") 
    include("common/timers.jl")
    include("common/get_active_nprocs.jl")
    include("common/print_results.jl")
    include("SP32/classes.jl")
    include("SP32/sp_data.jl")
    include("SP32/mpinpb.jl")
    include("SP32/setup_mpi.jl") 
    include("SP32/zone_setup.jl")
    include("SP32/make_set.jl") 
    include("SP32/define.jl") 
    include("SP32/initialize.jl") 
    include("SP32/set_constants.jl") 
    include("SP32/add.jl") 
    include("SP32/copy_faces.jl")  
    include("SP32/exchange_qbc.jl")  
    include("SP32/error.jl") 
    include("SP32/exact_rhs.jl") 
    include("SP32/exact_solution.jl") 
    include("SP32/lhsx.jl") 
    include("SP32/lhsy.jl") 
    include("SP32/lhsz.jl") 
    include("SP32/ninvr.jl") 
    include("SP32/pinvr.jl") 
    include("SP32/rhs.jl") 
    include("SP32/txinvr.jl") 
    include("SP32/tzetar.jl") 
    include("SP32/verify.jl") 
    include("SP32/x_solve.jl") 
    include("SP32/y_solve.jl") 
    include("SP32/z_solve.jl")
    include("SP32/adi.jl") 
    include("SP32/sp-node.jl") 
    include("SP32/sp-cluster.jl") 
    include("SP32/sp-driver.jl") 
end

module BT
    using FortranFiles
    using OffsetArrays
    using Parameters
    using Printf
    using MPI
    using StaticArrays
    using MPIClusterManagers
    using Distributed
 
    include("common/DFVariable.jl") 
    include("common/timers.jl")
    include("common/get_active_nprocs.jl")
    include("common/print_results.jl")
    include("BT32/classes.jl")
    include("BT32/bt_data.jl")
    include("BT32/set_constants.jl")
    include("BT32/mpinpb.jl")
    include("BT32/setup_mpi.jl")
    include("BT32/zone_setup.jl")
    include("BT32/make_set.jl")
    include("BT32/define.jl")
    include("BT32/initialize.jl")
    include("BT32/add.jl")
    include("BT32/btio.jl")
    include("BT32/copy_faces.jl")
    include("BT32/exchange_qbc.jl") 
    include("BT32/error.jl")
    include("BT32/exact_rhs.jl")
    include("BT32/exact_solution.jl")
    include("BT32/rhs.jl")
    include("BT32/solve_subs.jl")
    include("BT32/x_solve.jl")
    include("BT32/y_solve.jl")
    include("BT32/z_solve.jl")    
    include("BT32/verify.jl")
    include("BT32/adi.jl")
    include("BT32/bt-node.jl") 
    include("BT32/bt-cluster.jl") 
    include("BT32/bt-driver.jl") 

end

module LU
    using FortranFiles
    using OffsetArrays
    using Parameters
    using Printf
    using MPI
    using StaticArrays
    using MPIClusterManagers
    using Distributed
    using ProgressMeter

    include("common/DFVariable.jl") 
    include("common/timers.jl")
    include("common/get_active_nprocs.jl")
    include("common/print_results.jl")
    include("LU/classes.jl")
    include("LU/nodedim.jl")
    include("LU/neighbors.jl")
    include("LU/subdomain.jl")
    include("LU/setcoeff.jl")
#    include("LU/lu_data_vec.jl")
    include("LU/lu_data.jl")
    include("LU/init_comm.jl")
    include("LU/mpinpb.jl")
    include("LU/zone_setup.jl")
    include("LU/bcast_inputs.jl")
    include("LU/exchange_qbc.jl")
#    include("LU/blts_vec.jl")
    include("LU/blts.jl")
#    include("LU/buts_vec.jl")
    include("LU/buts.jl")
    include("LU/erhs.jl")
    include("LU/error.jl")
    include("LU/exact.jl")
    include("LU/exchange_1.jl")
    include("LU/exchange_3.jl")
    include("LU/exchange_4.jl")
    include("LU/exchange_5.jl")
    include("LU/exchange_6.jl")
    #include("LU/jacld_vec.jl")
    include("LU/jacld.jl")
#    include("LU/jacu_vec.jl")
    include("LU/jacu.jl")
    include("LU/l2norm.jl")
    include("LU/lu-node.jl")
    include("LU/lu-cluster.jl")
    include("LU/lu-driver.jl")
    include("LU/pintgr.jl")
    include("LU/proc_grid.jl")
    include("LU/rhs.jl")
    include("LU/setbv.jl")
    include("LU/setiv.jl")
#    include("LU/ssor_vec.jl")
    include("LU/ssor.jl")
    include("LU/verify.jl")

end

export SP, BT, LU

end # module NPBApps