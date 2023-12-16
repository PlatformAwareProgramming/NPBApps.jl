using FortranFiles
using OffsetArrays
using Parameters
using Printf


#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function read_input(class)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#      use lu_data
#      use mpinpb
#      use timing

#      implicit none

#      character class
#      integer IERROR, fstatus


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

         check_timer_flag( timeron )

         OPEN(unit = 3, file = "inputlu.data", status = "old",
               access = "sequential", form = "formatted", iostat = fstatus)
         if fstatus == 0

            println(stdout, "Reading from input file inputlu.data")

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
         else
            ipr = ipr_default
            inorm = inorm_default
            itmax = itmax_default
            dt = dt_default
            omega = omega_default
            tolrsd[1] = tolrsddef(1)
            tolrsd[2] = tolrsddef(2)
            tolrsd[3] = tolrsddef(3)
            tolrsd[4] = tolrsddef(4)
            tolrsd[5] = tolrsddef(5)
            nx0 = isiz01
            ny0 = isiz02
            nz0 = isiz03
         end

#---------------------------------------------------------------------
#   check problem size
#---------------------------------------------------------------------
         if ( nx0 < 4 ) ||(
               ny0 < 4 ) ||(
               nz0 < 4 )

            @printf(stdout, "     PROBLEM SIZE IS TOO SMALL - \n     SET EACH OF NX, NY AND NZ AT LEAST EQUAL TO 5\n", )
# 2001       format(5x,'PROBLEM SIZE IS TOO SMALL - ',  		                 /5x,'SET EACH OF NX, NY AND NZ AT LEAST EQUAL TO 5')
            MPI_ABORT( MPI_COMM_WORLD, MPI_ERR_OTHER, IERROR )

         end

         if ( nx0 > isiz01 ) ||(
               ny0 > isiz02 ) ||(
               nz0 > isiz03 )

            @printf(stdout, "     PROBLEM SIZE IS TOO LARGE - \n     NX, NY AND NZ SHOULD BE LESS THAN OR EQUAL TO \n     ISIZ01, ISIZ02 AND ISIZ03 RESPECTIVELY\n", )
# 2002       format(5x,'PROBLEM SIZE IS TOO LARGE - ',  		                 /5x,'NX, NY AND NZ SHOULD BE LESS THAN OR EQUAL TO ',  		                 /5x,'ISIZ01, ISIZ02 AND ISIZ03 RESPECTIVELY')
            MPI_ABORT( MPI_COMM_WORLD, MPI_ERR_OTHER, IERROR )

         end

         set_class(class)

         @printf(stdout, " Size: %4ix%4ix%4i  (class %s)\n", nx0, ny0, nz0, class)
         @printf(stdout, " Iterations: %4i\n", itmax)

         @printf(stdout, " Total number of processes: %6i\n", total_nodes)
         if (total_nodes != no_nodes) @printf(stdout, " WARNING: Number of processes is not in a form of (n1*n2, n1/n2 <= 2).\n Number of active processes: %6i\n", no_nodes) end
         println(stdout, )

# 1000 format(//, ' NAS Parallel Benchmarks 3.4 -- LU Benchmark',/)
# 1001    format(' Size: ', i4, 'x', i4, 'x', i4, '  (class ', a, ')')
# 1002    format(' Iterations: ', i4)
# 1003    format(' Total number of processes: ', i6)
# 1004    format(' WARNING: Number of processes is not in a form of',  		                ' (n1*n2, n1/n2 <= 2).'/  		                ' Number of active processes: ', i6)


      end

      bcast_inputs

      return nothing
      end


