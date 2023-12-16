module NPB_MZ_Multicluster

using FortranFiles
using OffsetArrays
using Parameters
using Printf
using MPI

include("common/timers.jl")
include("common/get_active_nprocs.jl")
include("common/print_results.jl")
include("SP/npbparams.jl")
include("SP/sp_data.jl")
# include("SP/mpinpb.jl") 
include("SP/setup_mpi.jl") 
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
include("SP/sp.jl") 



export SP

end # module NPB_MZ_Multicluster