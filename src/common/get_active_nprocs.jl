

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# return the largest np1 and np2 such that np1 * np2 <= nprocs
# pkind = 1, np1 = np2                 (square number)
#         2, np1/2 <= np2 <= np1
#         3, np1 = np2 or np1 = np2*2  (power of 2)
# other outputs:
#     npa = np1 * np2 (active number of processes)
#     nprocs   - total number of processes
#     rank     - rank of this process
#     comm_out - MPI communicator
#     active   - .true. if this process is active; .false. otherwise
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function get_active_nprocs(comm_in, pkind)

# nprocs and rank in comm_in
      nprocs = MPI.Comm_size(comm_in) # ios = ?
      RANK = MPI.Comm_rank(comm_in)   # ios = ?
      
      if pkind <= 1
# square number of processes (add small number to allow for roundoff)
         np2 = trunc(Int, sqrt(float(nprocs) + 1.0e-3))
         np1 = np2
      else
# power-of-two processes
         np1 = trunc(Int, log(float(nprocs) + 1.0e-3) / log(2.0e0))
         np2 = np1 / 2
         np1 = np1 - np2
         np1 = 2^np1
         np2 = 2^np2
      end
      npa = np1 * np2

# for option 2, go further to get the best (np1 * np2) proc grid
      if pkind == 2 && npa < nprocs
         np1w = trunc(Int, sqrt(float(nprocs) + 1.0e-3))
         np2w = trunc(Int, sqrt(float(nprocs*2) + 1.0e-3))
         for n = np1w:np2w
            npaw = nprocs / n * n
            if (n == 1 && nprocs == 3) npaw = 2 end
            if npaw > npa
               npa = npaw
               np1 = npa / n
               if np1 < n
                  np2 = np1
                  np1 = n
               else
                  np2 = n
               end
            end
         end
      end

# all good if calculated is the same as requested
      comm_out = comm_in
      active = true
      if (nprocs == npa) 
         return Int(np1), Int(np2), Int(npa), nprocs, RANK, comm_out, active
      end

# npa < nprocs, need to check if a strict NPROCS enforcement is required
      if RANK == 0
         ios = haskey(ENV, "NPB_NPROCS_STRICT") ? 0 : 1
         if ios == 0 #&& ic > 0
            val = ENV["NPB_NPROCS_STRICT"]
            if val == "0" || val[1] == '-'
               active = false
            elseif val == "off" || val == "OFF" ||
                  val[1] == 'n' || val[1] == 'N' ||
                  val[1] == 'f' || val[1] == 'F'
               active = false
            end
         end
      end

      active = MPI.Bcast(true, 0, comm_in)

# abort if a strict NPROCS enforcement is required
      if active
         if RANK == 0
            println(" *** ERROR determining processor topology for $nprocs processes")
            if pkind <= 1
               t = "square"
            elseif pkind == 2
               t = "grid (nx*ny, nx/2<=ny<=nx)"
            else
               t = "power-of-two"
            end
            println("     Expecting a $t number of processes (such as $npa)")
         end
         ios = MPI.Abort(comm_in, MPI.MPI_ERR_OTHER) # ios = ?
         exit(0)
      end

# mark excess ranks as inactive
# split communicator based on rank value
      if RANK >= npa
         active = false
         ic = 1
      else
         active = true
         ic = 0
      end

      comm_out = MPI.Comm_split(comm_in, ic, RANK)  # ios = ?
      
      return Int(np1), Int(np2), Int(npa), nprocs, RANK, comm_out, activev

end
