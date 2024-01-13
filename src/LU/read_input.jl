


#---------------------------------------------------------------------
#---------------------------------------------------------------------

 function read_input(params_file)

#---------------------------------------------------------------------
#    only root reads the input file
#    if input file does not exist, it uses defaults
#       ipr = 1 for detailed progress output
#       inorm = how often the norm is printed (once every inorm iterations)
#       itmax = number of pseudo time steps
#       dt = time step
#       omega 1 over-relaxation factor for SSOR
#       tolrsd = steady state residual tolerance levels
#       nx, ny, nz = number of grid points in x, y, z directions
#---------------------------------------------------------------------

      if id == root

         @printf(stdout, "\n\n NAS Parallel Benchmarks 3.4 -- LU Benchmark\n\n", )

         timeron = check_timer_flag()

         #OPEN(unit = 3, file = "inputlu.data", status = "old",      access = "sequential", form = "formatted", iostat = fstatus)

         fstatus = isfile(params_file) ? 0 : 1      
         #=if fstatus == 0
            println(stdout, "Reading from input file inputlu.data")
            f = open("inputlu.data", "r")

            READ(3, )
            READ(3, )
            READ(3, ipr, inorm)
            READ(3, )
            READ(3, )
            READ(3, itmax)
            READ(3, )
            READ(3, )
            READ(3, dt)
            READ(3, )
            READ(3, )
            READ(3, omega)
            READ(3, )
            READ(3, )
            READ(3, tolrsd[1], tolrsd[2], tolrsd[3], tolrsd[4], tolrsd[5])
            READ(3, )
            READ(3, )
            READ(3, nx0, ny0, nz0)
            CLOSE(3)
         else=#
         #end

         # class = set_class(itmax, isiz01, isiz02, isiz03)

      end

      #bcast_inputs()

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
      
      return nx0, nx1, nx2, itmax, inorm, dt
end


