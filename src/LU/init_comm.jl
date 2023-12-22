#---------------------------------------------------------------------
#
#   initialize MPI and establish rank and size
#
# This is a module in the MPI implementation of LUSSOR
# pseudo application from the NAS Parallel Benchmarks. 
#
#---------------------------------------------------------------------

 function init_comm()

#---------------------------------------------------------------------
#    initialize MPI communication
#---------------------------------------------------------------------
      MPI.Init()

#---------------------------------------------------------------------
#     get a process grid that requires a (nx*ny) number of procs.
#     excess ranks are marked as inactive.
#---------------------------------------------------------------------

      global xdim, ydim, no_nodes, total_nodes, node, comm_solve, active = get_active_nprocs(MPI.COMM_WORLD, 2)

      if (!active) return end

#---------------------------------------------------------------------
#   establish the global rank of this process and the group size
#---------------------------------------------------------------------
      global id = node
      global num = no_nodes
      global root = 0

      global ndim   = nodedim(num)

 #=     if !convertdouble
         dp_type = MPI_DOUBLE_PRECISION
      else
         dp_type = MPI_REAL
      end
=#

      return nothing
end
