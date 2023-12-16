

function print_results(name, class, n1, n2, n3, niter,
                     nprocs_active, nprocs_total,
                     t, mops, optype, verified, npbversion,
                     compiletime, cs1, cs2, cs3, cs4, cs5, cs6, cs7)

         @printf(stdout, "\n\n %2s Benchmark Completed.\n", name)

         @printf(stdout, " Class           =             %12s\n", class)

         if (n2 == 0) && (n3 == 0)
            if name[1:2] == "EP"
               @printf(SIZE, "%15.0F\n", 2.0e0^n1)
               j = 15
               if (SIZE[j:j] == ".") j = j - 1 end
               @printf(stdout, " Size            =          %15s\n", SIZE[1:j])
            else
               @printf(stdout, " Size            =             %12i\n", n1)
            end
         else
            @printf(stdout, " Size            =           %4ix%4ix%4i\n", n1, n2, n3)
         end

         @printf(stdout, " Iterations      =             %12i\n", niter)
         @printf(stdout, " Time in seconds =             %12.2F\n", t)
         @printf(stdout, " Total processes =             %12i\n", nprocs_total)

         if (nprocs_active != 0) 
            @printf(stdout, " Active processes=             %12i\n", nprocs_active) 
         end

         @printf(stdout, " Mop/s total     =             %12.2F\n", mops)
         @printf(stdout, " Mop/s/process   =             %12.2F\n", mops/float( nprocs_total ))
         @printf(stdout, " Operation type  = %24s\n", optype)

         if verified
            @printf(stdout, " Verification    =             %s\n", "  SUCCESSFUL")
         else
            @printf(stdout, " Verification    =             %s\n", "UNSUCCESSFUL")
         end

         @printf(stdout, " Version         =             %12s\n", npbversion)
#         @printf(stdout, " Compile date    =             %12s\n", compiletime)
#         @printf(stdout, "\n Compile options:\n    MPIFC        = %s\n", cs1)
#         @printf(stdout, "    FLINK        = %s\n", cs2)
#         @printf(stdout, "    FMPI_LIB     = %s\n", cs3)
#         @printf(stdout, "    FMPI_INC     = %s\n", cs4)
#         @printf(stdout, "    FFLAGS       = %s\n", cs5)
#         @printf(stdout, "    FLINKFLAGS   = %s\n", cs6)
#         @printf(stdout, "    RAND         = %s\n", cs7)
         @printf(stdout, "\n\n Please send feedbacks and/or the results of this run to:\n\n NPB Development Team \n Internet: npb@nas.nasa.gov\n\n\n", )

         return nothing
end

