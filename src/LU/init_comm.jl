#---------------------------------------------------------------------
#
#   initialize MPI and establish rank and size
#
# This is a module in the MPI implementation of LUSSOR
# pseudo application from the NAS Parallel Benchmarks. 
#
#---------------------------------------------------------------------
 const root = 0

 function init_comm()

#---------------------------------------------------------------------
#    initialize MPI communication
#---------------------------------------------------------------------
      MPI.Init()

#---------------------------------------------------------------------
#     get a process grid that requires a (nx*ny) number of procs.
#     excess ranks are marked as inactive.
#---------------------------------------------------------------------

      xdim, ydim, no_nodes, total_nodes, node, comm_solve, active = get_active_nprocs(MPI.COMM_WORLD, 2)

      if (!active) return end

#---------------------------------------------------------------------
#   establish the global rank of this process and the group size
#---------------------------------------------------------------------
      id = node
      num = no_nodes

      ndim   = nodedim(num)

 #=     if !convertdouble
         dp_type = MPI_DOUBLE_PRECISION
      else
         dp_type = MPI_REAL
      end
=#

      return xdim, ydim, no_nodes, total_nodes, node, comm_solve, active, id, num, ndim
end
