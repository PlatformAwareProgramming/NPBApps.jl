


#---------------------------------------------------------------------
#---------------------------------------------------------------------

 function read_input()

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

      global ipr = ipr_default
      global inorm = inorm_default
      global itmax = itmax_default
      global dt = dt_default
      global omega = omega_default
      global tolrsd = tolrsddef
      global nx0 = isiz01
      global ny0 = isiz02
      global nz0 = isiz03
      global timeron = 0
      class = "U"

      if id == root

         @printf(stdout, "\n\n NAS Parallel Benchmarks 3.4 -- LU Benchmark\n\n", )

         timeron = check_timer_flag()

         #OPEN(unit = 3, file = "inputlu.data", status = "old",      access = "sequential", form = "formatted", iostat = fstatus)

         fstatus = isfile("inputlu.data") ? 0 : 1      
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

#---------------------------------------------------------------------
#   check problem size
#---------------------------------------------------------------------
         if (nx0 < 4) || (ny0 < 4 ) || (nz0 < 4)

            @printf(stdout, "     PROBLEM SIZE IS TOO SMALL - \n     SET EACH OF NX, NY AND NZ AT LEAST EQUAL TO 5\n", )
            MPI.Abort(MPI.COMM_WORLD, MPI.MPI_ERR_OTHER)

         end

         if (nx0 > isiz01) || (ny0 > isiz02) || (nz0 > isiz03)

            @printf(stdout, "     PROBLEM SIZE IS TOO LARGE - \n     NX, NY AND NZ SHOULD BE LESS THAN OR EQUAL TO \n     ISIZ01, ISIZ02 AND ISIZ03 RESPECTIVELY\n", )
            MPI.Abort(MPI.COMM_WORLD, MPI.MPI_ERR_OTHER)

         end

         class = set_class()

         @printf(stdout, " Size: %4ix%4ix%4i  (class %s)\n", nx0, ny0, nz0, class)
         @printf(stdout, " Iterations: %4i\n", itmax)

         @printf(stdout, " Total number of processes: %6i\n", total_nodes)
         if (total_nodes != no_nodes) @printf(stdout, " WARNING: Number of processes is not in a form of (n1*n2, n1/n2 <= 2).\n Number of active processes: %6i\n", no_nodes) end
         println(stdout, )

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
      
      return class
end


