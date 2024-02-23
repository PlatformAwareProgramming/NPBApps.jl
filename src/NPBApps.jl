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
    include("SP/classes.jl")
    include("SP/sp_data.jl")
    include("SP/mpinpb.jl")
    include("SP/setup_mpi.jl") 
    include("SP/zone_setup.jl")
    include("SP/make_set.jl") 
    include("SP/define.jl") 
    include("SP/initialize.jl") 
    include("SP/set_constants.jl") 
    include("SP/add.jl") 
    include("SP/copy_faces.jl") 
    include("SP/error.jl") 
    include("SP/exact_rhs.jl") 
    include("SP/exact_solution.jl") 
    include("SP/lhsx.jl") 
    include("SP/lhsy.jl") 
    include("SP/lhsz.jl") 
    include("SP/ninvr.jl") 
    include("SP/pinvr.jl") 
    include("SP/rhs.jl") 
    include("SP/txinvr.jl") 
    include("SP/tzetar.jl") 
    include("SP/verify.jl") 
    include("SP/x_solve.jl") 
    include("SP/y_solve.jl") 
    include("SP/z_solve.jl")
    include("SP/adi.jl") 
    include("SP/sp-node.jl") 
    include("SP/sp-cluster.jl") 
    include("SP/sp-driver.jl") 
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
    include("BT/classes.jl")
    include("BT/bt_data.jl")
    include("BT/set_constants.jl")
    include("BT/mpinpb.jl")
    include("BT/setup_mpi.jl")
    include("BT/zone_setup.jl")
    include("BT/make_set.jl")
    include("BT/define.jl")
    include("BT/initialize.jl")
    include("BT/add.jl")
    include("BT/btio.jl")
    include("BT/copy_faces.jl")
    include("BT/error.jl")
    include("BT/exact_rhs.jl")
    include("BT/exact_solution.jl")
    include("BT/rhs.jl")
    include("BT/solve_subs.jl")
    include("BT/x_solve.jl")
    include("BT/y_solve.jl")
    include("BT/z_solve.jl")    
    include("BT/verify.jl")
    include("BT/adi.jl")
    include("BT/bt-node.jl") 
    include("BT/bt-cluster.jl") 
    include("BT/bt-driver.jl") 

end

module LU
    using FortranFiles
    using OffsetArrays
    using Parameters
    using Printf
    using MPI
    using StaticArrays

    include("common/timers.jl")
    include("common/get_active_nprocs.jl")
    include("common/print_results.jl")
    include("LU/classes.jl")
    include("LU/npbparams.jl")
    include("LU/nodedim.jl")
    include("LU/neighbors.jl")
    include("LU/subdomain.jl")
    include("LU/setcoeff.jl")
#    include("LU/lu_data_vec.jl")
    include("LU/lu_data.jl")
    include("LU/init_comm.jl")
    include("LU/read_input.jl")
    include("LU/bcast_inputs.jl")
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
    include("LU/jacld_vec.jl")
    include("LU/jacld.jl")
#    include("LU/jacu_vec.jl")
    include("LU/jacu.jl")
    include("LU/l2norm.jl")
    include("LU/lu.jl")
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