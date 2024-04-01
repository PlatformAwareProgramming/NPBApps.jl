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

const t_recs = "total", "rhs", "blts", "buts", "#jacld", "#jacu",
            "exch", "lcomm", "ucomm", "rcomm",
            " totcomp", " totcomm"

function go(class::CLASS; timers=false)
         
      itmax = bt_class[class].itmax
      inorm = bt_class[class].inorm
      dt    = bt_class[class].dt
      isiz01 = bt_class[class].isiz01
      isiz02 = bt_class[class].isiz02
      isiz03 = bt_class[class].isiz03
   
      go(isiz01, isiz02, isiz03, itmax, inorm, dt; timers=timers)
   
end
   
function go(isiz01, isiz02, isiz03, itmax, inorm, dt; timers=false)

      xdim, ydim, no_nodes, total_nodes, node, comm_solve, active, id, num, ndim = init_comm()

      perform(xdim, ydim, no_nodes, total_nodes, node, comm_solve, active, id, num, ndim, isiz01, isiz02, isiz03, itmax, inorm, dt, timers)

end

function perform(xdim, ydim, no_nodes, total_nodes, node, comm_solve, active, id, num, ndim, isiz01, isiz02, isiz03, itmax, inorm, dt, timeron)   


      class = set_class(itmax, isiz01, isiz02, isiz03)

#---------------------------------------------------------------------
#   initialize communications
#---------------------------------------------------------------------
      if (!active) @goto L999 end

#---------------------------------------------------------------------
#   read input data
#---------------------------------------------------------------------
      
      omega = omega_default
      tolrsd = tolrsddef
      nx0 = isiz01
      ny0 = isiz02
      nz0 = isiz03

      if id == root

         @printf(stdout, "\n\n NAS Parallel Benchmarks 3.4 -- LU Benchmark\n\n", )

#---------------------------------------------------------------------
#   check problem size
#---------------------------------------------------------------------
         if (nx0 < 4) || (ny0 < 4 ) || (nz0 < 4)

            @printf(stdout, "     PROBLEM SIZE IS TOO SMALL - \n     SET EACH OF NX, NY AND NZ AT LEAST EQUAL TO 5\n", )
            MPI.Abort(MPI.COMM_WORLD, MPI.MPI_ERR_OTHER)

         end

         #=if (nx0 > isiz01) || (ny0 > isiz02) || (nz0 > isiz03)

            @printf(stdout, "     PROBLEM SIZE IS TOO LARGE - \n     NX, NY AND NZ SHOULD BE LESS THAN OR EQUAL TO \n     ISIZ01, ISIZ02 AND ISIZ03 RESPECTIVELY\n", )
            MPI.Abort(MPI.COMM_WORLD, MPI.MPI_ERR_OTHER)

         end=#

         @printf(stdout, " Size: %4ix%4ix%4i  (class %s)\n", nx0, ny0, nz0, class)
         @printf(stdout, " Iterations: %4i\n", itmax)

         @printf(stdout, " Total number of processes: %6i\n", total_nodes)
         if (total_nodes != no_nodes) @printf(stdout, " WARNING: Number of processes is not in a form of (n1*n2, n1/n2 <= 2).\n Number of active processes: %6i\n", no_nodes) end
         println(stdout, )

      end

      for i = 1:t_last
         timer_clear(i)
      end

#---------------------------------------------------------------------
#   set up processor grid
#---------------------------------------------------------------------
      isiz1, isiz2, isiz3, row, col = proc_grid(id, nx0, ny0, nz0, no_nodes, num, ndim)

#---------------------------------------------------------------------
#   allocate space
#---------------------------------------------------------------------
      rsdnm, errnm, u, rsd, frct, flux, a, b, c, d, phi1, phi2, buf, buf1 = alloc_space(isiz1, isiz2, isiz3)

#---------------------------------------------------------------------
#   determine the neighbors
#---------------------------------------------------------------------
     south, east, north, west = neighbors(row, col, xdim, ydim)

#---------------------------------------------------------------------
#   set up sub-domain sizes
#---------------------------------------------------------------------
      nx, ny, nz, ipt, jpt, ist, jst, iend, jend = subdomain(row, col, nx0, ny0, nz0, isiz1, isiz2, isiz3, north, south, west, east, xdim, ydim)

#---------------------------------------------------------------------
#   set up coefficients
#---------------------------------------------------------------------
      dxi, deta, dzeta, tx1, tx2, tx3, ty1, ty2, ty3, tz1, tz2, tz3, ii2, ji2, ki2 = setcoeff(nx0, ny0, nz0)

#---------------------------------------------------------------------
#   set the boundary values for dependent variables
#---------------------------------------------------------------------
      setbv(u, nx0, ny0, nz0, nx, ny, nz, ipt, jpt, north, south, west, east)

#---------------------------------------------------------------------
#   set the initial values for dependent variables
#---------------------------------------------------------------------
      setiv(u, nx0, ny0, nz0, nx, ny, nz, ipt, jpt)

#---------------------------------------------------------------------
#   compute the forcing term based on prescribed exact solution
#---------------------------------------------------------------------
      erhs(frct, rsd, flux, buf, buf1, nx0, ny0, nx, ny, nz, ipt, jpt, ist, jst, iend, jend, north, south, west, east, 
           dssp, tx1, tx2, tx3, ty1, ty2, ty3, tz1, tz2, tz3, 
           comm_solve)

#---------------------------------------------------------------------
#   perform one SSOR iteration to touch all data and program pages 
#---------------------------------------------------------------------
      ssor(1,
            comm_solve, 
            id,
            rsdnm,
            u,
            rsd,
            frct,
            flux,
            a,
            b,
            c,
            d,
            buf,
            buf1,
            south,
            east,
            north,
            west,
            isiz1,
            isiz2,
            isiz3,
            inorm,
            itmax,
            dt,
            omega,
            tolrsd,
            nx0,
            ny0,
            nz0,
            timeron,
            tx1,
            tx2,
            tx3,
            ty1,
            ty2,
            ty3,
            tz1,
            tz2,
            tz3,
            nx,
            ipt,
            ny,
            jpt,
            nz,
            ist,
            iend,
            jst,
            jend,
      )
#---------------------------------------------------------------------
#   reset the boundary and initial values
#---------------------------------------------------------------------
      setbv(u, nx0, ny0, nz0, nx, ny, nz, ipt, jpt, north, south, west, east)
      setiv(u, nx0, ny0, nz0, nx, ny, nz, ipt, jpt)

#---------------------------------------------------------------------
#   perform the SSOR iterations
#---------------------------------------------------------------------
      maxtime = ssor(
                  itmax,
                  comm_solve, 
                  id,
                  rsdnm,
                  u,
                  rsd,
                  frct,
                  flux,
                  a,
                  b,
                  c,
                  d,
                  buf,
                  buf1,
                  south,
                  east,
                  north,
                  west,
                  isiz1,
                  isiz2,
                  isiz3,
                  inorm,
                  itmax,
                  dt,
                  omega,
                  tolrsd,
                  nx0,
                  ny0,
                  nz0,
                  timeron,
                  tx1,
                  tx2,
                  tx3,
                  ty1,
                  ty2,
                  ty3,
                  tz1,
                  tz2,
                  tz3,
                  nx,
                  ipt,
                  ny,
                  jpt,
                  nz,
                  ist,
                  iend,
                  jst,
                  jend,
            )

#---------------------------------------------------------------------
#   compute the solution error
#---------------------------------------------------------------------
      ERROR(u, errnm, nx0, ny0, nz0, nz, jst, jpt, ist, ipt, iend, jend, comm_solve)

#---------------------------------------------------------------------
#   compute the surface integral
#---------------------------------------------------------------------
      frc = pintgr(u, nx, ny, nz, isiz02, isiz03, isiz2, isiz3, ipt, jpt, north, south, west, east, phi1, phi2,  dxi, deta, dzeta, ii2, ji2, ki2, comm_solve)
      
#---------------------------------------------------------------------
#   verification test
#---------------------------------------------------------------------
      if id == 0
         verified = verify( rsdnm, errnm, frc, class, dt)
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

      if timeron
         t1 = Array{Float64}(undef, t_last+2)
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
      end
      @label L999
      MPI.Finalize()

     return nothing
end


