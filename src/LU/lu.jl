#-------------------------------------------------------------------------!
#                                                                         !
#        N  A  S     P A R A L L E L     B E N C H M A R K S  3.4         !
#                                                                         !
#                                   L U                                   !
#                                                                         !
#-------------------------------------------------------------------------!
#                                                                         !
#    This benchmark is part of the NAS Parallel Benchmark 3.4 suite.      !
#    It is described in NAS Technical Reports 95-020 and 02-007           !
#                                                                         !
#    Permission to use, copy, distribute and modify this software         !
#    for any purpose with or without fee is hereby granted.  We           !
#    request, however, that all derived work reference the NAS            !
#    Parallel Benchmarks 3.4. This software is provided "as is"           !
#    without express or implied warranty.                                 !
#                                                                         !
#    Information on NPB 3.4, including the technical report, the          !
#    original specifications, source code, results and information        !
#    on how to submit new results, is available at:                       !
#                                                                         !
#           http://www.nas.nasa.gov/Software/NPB/                         !
#                                                                         !
#    Send comments or suggestions to  npb@nas.nasa.gov                    !
#                                                                         !
#          NAS Parallel Benchmarks Group                                  !
#          NASA Ames Research Center                                      !
#          Mail Stop: T27A-1                                              !
#          Moffett Field, CA   94035-1000                                 !
#                                                                         !
#          E-mail:  npb@nas.nasa.gov                                      !
#          Fax:     (650) 604-3957                                        !
#                                                                         !
#-------------------------------------------------------------------------!

#---------------------------------------------------------------------
#
# Authors: S. Weeratunga
#          V. Venkatakrishnan
#          E. Barszcz
#          M. Yarrow
#
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#
#   driver for the performance evaluation of the solver for
#   five coupled parabolic/elliptic partial differential equations.
#
#---------------------------------------------------------------------

#---------------------------------------------------------------------
 function go()
#---------------------------------------------------------------------

      t_recs = "total", "rhs", "blts", "buts", "#jacld", "#jacu",
            "exch", "lcomm", "ucomm", "rcomm",
            " totcomp", " totcomm"

#---------------------------------------------------------------------
#   initialize communications
#---------------------------------------------------------------------
      init_comm()
      if (!active) @goto L999 end

#---------------------------------------------------------------------
#   read input data
#---------------------------------------------------------------------
      class = read_input()

      for i = 1:t_last
         timer_clear(i)
      end

#---------------------------------------------------------------------
#   set up processor grid
#---------------------------------------------------------------------
      proc_grid()

#---------------------------------------------------------------------
#   allocate space
#---------------------------------------------------------------------
      alloc_space()

#---------------------------------------------------------------------
#   determine the neighbors
#---------------------------------------------------------------------
      neighbors()

#---------------------------------------------------------------------
#   set up sub-domain sizes
#---------------------------------------------------------------------
      subdomain()

#---------------------------------------------------------------------
#   set up coefficients
#---------------------------------------------------------------------
      setcoeff()

#---------------------------------------------------------------------
#   set the boundary values for dependent variables
#---------------------------------------------------------------------
      setbv()

#---------------------------------------------------------------------
#   set the initial values for dependent variables
#---------------------------------------------------------------------
      setiv()

#---------------------------------------------------------------------
#   compute the forcing term based on prescribed exact solution
#---------------------------------------------------------------------
      erhs()

#---------------------------------------------------------------------
#   perform one SSOR iteration to touch all data and program pages 
#---------------------------------------------------------------------
      ssor(1)

#---------------------------------------------------------------------
#   reset the boundary and initial values
#---------------------------------------------------------------------
      setbv()
      setiv()

#---------------------------------------------------------------------
#   perform the SSOR iterations
#---------------------------------------------------------------------
      ssor(itmax)

#---------------------------------------------------------------------
#   compute the solution error
#---------------------------------------------------------------------
      ERROR()

#---------------------------------------------------------------------
#   compute the surface integral
#---------------------------------------------------------------------
      frc = pintgr()
      
#---------------------------------------------------------------------
#   verification test
#---------------------------------------------------------------------
      if id == 0
         verified = verify( rsdnm, errnm, frc, class)
         mflops = 1.0e-6*float(itmax)*(1984.77*float( nx0 )*
              float( ny0 )*
              float( nz0 )-
              10923.3*(float( nx0+ny0+nz0 )/3.)^2+
              27770.9* float( nx0+ny0+nz0 )/3.0-
              144010.)/
               maxtime

         print_results("LU", class, nx0,
           ny0, nz0, itmax, no_nodes, total_nodes,
           maxtime, mflops, "          floating point", verified,
           npbversion)

      end

      if (!timeron) @goto L999 end

      for i = 1:t_last
         t1[i] = timer_read(i)
      end
      t1[t_rhs] = t1[t_rhs] - t1[t_exch]
      t1[t_last+2] = t1[t_lcomm]+t1[t_ucomm]+t1[t_rcomm]+t1[t_exch]
      t1[t_last+1] = t1[t_total] - t1[t_last+2]

      tsum = MPI.Reduce(t1, MPI.SUM, 0, comm_solve)
      tming = MPI.Reduce(t1, MPI.MIN, 0, comm_solve)
      tmaxg = MPI.Reduce(t1, MPI.MAX, 0, comm_solve)

      if id == 0
         @printf(stdout, " nprocs =%6i           minimum     maximum     average\n", no_nodes)
         for i = 1:t_last+2
            if t_recs[i][1:1] != "#"
               tsum[i] = tsum[i] / no_nodes
               @printf(stdout, " timer %2i(%8s) :  %10.4F  %10.4F  %10.4F\n", i, t_recs[i], tming[i], tmaxg[i], tsum[i])
            end
         end
      end

      @label L999
      MPI.Finalize()

     return nothing
end


