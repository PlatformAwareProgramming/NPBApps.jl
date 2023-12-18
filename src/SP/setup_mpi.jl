

#---------------------------------------------------------------------
# set up MPI stuff
#---------------------------------------------------------------------

function setup_mpi()

      global DEFAULT_TAG = 0

      MPI.Init()

#---------------------------------------------------------------------
#     get a process grid that requires a square number of procs.
#     excess ranks are marked as inactive.
#---------------------------------------------------------------------
      global _, maxcells, no_nodes, total_nodes, node, comm_setup, active = get_active_nprocs(MPI.COMM_WORLD, 1)

      if (!active) return end

      global comm_solve = MPI.Comm_dup(comm_setup)
      global comm_rhs = MPI.Comm_dup(comm_setup)

#---------------------------------------------------------------------
#     let node 0 be the root for the group (there is only one)
#---------------------------------------------------------------------
      global root = 0

      return nothing
end

