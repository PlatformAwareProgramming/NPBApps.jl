using FortranFiles
using OffsetArrays
using Parameters
using Printf


#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function init_comm()

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#
#   initialize MPI and establish rank and size
#
# This is a module in the MPI implementation of LUSSOR
# pseudo application from the NAS Parallel Benchmarks. 
#
#---------------------------------------------------------------------

#      use lu_data
#      use mpinpb

#      implicit none

#      integer nodedim
#      integer IERROR


#---------------------------------------------------------------------
#    initialize MPI communication
#---------------------------------------------------------------------
      MPI_INIT( IERROR )

#---------------------------------------------------------------------
#     get a process grid that requires a (nx*ny) number of procs.
#     excess ranks are marked as inactive.
#---------------------------------------------------------------------
      get_active_nprocs(2, xdim, ydim, no_nodes,
                             total_nodes, node, comm_solve, active)

      if (!active) return end

#---------------------------------------------------------------------
#   establish the global rank of this process and the group size
#---------------------------------------------------------------------
      id = node
      num = no_nodes
      root = 0

      ndim   = nodedim(num)

      if !convertdouble
         dp_type = MPI_DOUBLE_PRECISION
      else
         dp_type = MPI_REAL
      end


      return nothing
      end
