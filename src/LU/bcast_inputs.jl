#---------------------------------------------------------------------
#---------------------------------------------------------------------

 function bcast_inputs()

#---------------------------------------------------------------------
#   root broadcasts the data
#   The data isn't contiguous or of the same type, so it's not
#   clear how to send it in the "MPI" way. 
#   We could pack the info into a buffer or we could create
#   an obscene datatype to handle it all at once. Since we only
#   broadcast the data once, just use a separate broadcast for
#   each piece. 
#---------------------------------------------------------------------
      ipr = MPI.bcast(ipr, comm_solve; root = root)
      inorm = MPI.bcast(inorm, comm_solve; root = root)
      itmax = MPI.bcast(itmax, comm_solve; root = root)
      dt = MPI.bcast(dt, comm_solve; root = root)
      omega = MPI.bcast(omega, comm_solve; root = root)
      tolrsd = MPI.bcast(tolrsd, comm_solve; root = root)
      nx0 = MPI.bcast(nx0, comm_solve; root = root)
      ny0 = MPI.bcast(ny0, comm_solve; root = root)
      nz0 = MPI.bcast(nz0, comm_solve; root = root)
      timeron = MPI.bcast(timeron, comm_solve; root = root)

      return nothing
end




