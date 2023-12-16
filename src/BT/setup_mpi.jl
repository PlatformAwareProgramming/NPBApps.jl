using FortranFiles
using OffsetArrays
using Parameters
using Printf


#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function setup_mpi()

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# set up MPI stuff
#---------------------------------------------------------------------

#      use bt_data
#      use mpinpb

#      implicit none

#      integer ERROR, nc

      mpi_init(ERROR)

      if !convertdouble
         dp_type = MPI_DOUBLE_PRECISION
      else
         dp_type = MPI_REAL
      end

#---------------------------------------------------------------------
#     get a process grid that requires a square number of procs.
#     excess ranks are marked as inactive.
#---------------------------------------------------------------------
      get_active_nprocs(1, nc, maxcells, no_nodes,
                             total_nodes, node, comm_setup, active)

      if (!active) return end

      mpi_comm_dup(comm_setup, comm_solve, ERROR)
      mpi_comm_dup(comm_setup, comm_rhs, ERROR)

#---------------------------------------------------------------------
#     let node 0 be the root for the group (there is only one)
#---------------------------------------------------------------------
      root = 0

      return nothing
      end

