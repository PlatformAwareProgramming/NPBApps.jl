#---------------------------------------------------------------------
# set up MPI stuff
#---------------------------------------------------------------------

const DEFAULT_TAG = 0

#---------------------------------------------------------------------
#     let node 0 be the root for the group (there is only one)
#---------------------------------------------------------------------
const root  = 0

function setup_mpi()

      
      MPI.Init()

#---------------------------------------------------------------------
#     get a process grid that requires a square number of procs.
#     excess ranks are marked as inactive.
#---------------------------------------------------------------------
      
      _, maxcells, no_nodes, total_nodes, node, comm_setup, active = get_active_nprocs(MPI.COMM_WORLD, 1)

      if (!active) return end

      comm_solve = MPI.Comm_dup(comm_setup)
      comm_rhs = MPI.Comm_dup(comm_setup)

      #return root, comm_solve, comm_rhs
      return maxcells, no_nodes, total_nodes, node, comm_setup, active, comm_solve, comm_rhs
end

