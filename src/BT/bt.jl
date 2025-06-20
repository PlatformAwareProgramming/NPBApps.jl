#-------------------------------------------------------------------------!
#                                                                         !
#        N  A  S     P A R A L L E L     B E N C H M A R K S  3.4         !
#                                                                         !
#                                   B T                                   !
#                                                                         !
#-------------------------------------------------------------------------!
#                                                                         !
#    This benchmark is part of the NAS Parallel Benchmark 3.4 suite.      !
#    It is described in NAS Technical Reports 95-020 and 02-007.          !
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
# Authors: R. F. Van der Wijngaart
#          T. Harris
#          M. Yarrow
#
#---------------------------------------------------------------------


function go(class::CLASS; timers = false, np=1, exeflags=``, mpiflags=``, dir="", threadlevel=:multiple, master_tcp_interface="")

   problem_size = bt_class[class].problem_size
   
   niter = bt_class[class].niter
   dt    = bt_class[class].dt

   grid_points = zeros(Integer, 3)
   grid_points[1] = problem_size
   grid_points[2] = problem_size
   grid_points[3] = problem_size

   workers = addprocs(MPIWorkerManager(np); threadlevel=threadlevel, exeflags=exeflags, mpiflags=mpiflags, dir=dir, master_tcp_interface=master_tcp_interface)

   try
      @everywhere workers @eval using NPBApps
      @everywhere workers BT.go($grid_points, $niter, $dt; timers = $timers)
   finally
      rmprocs(workers)
   end

end


function go(grid_points, niter, dt; timers = false)

   maxcells, no_nodes, total_nodes, node, comm_setup, active, comm_solve, comm_rhs = setup_mpi()

   perform(maxcells, no_nodes, total_nodes, node, comm_setup, active, comm_solve, comm_rhs, grid_points, niter, dt; timeron = timers)

end

function perform(maxcells, no_nodes, total_nodes, node, comm_setup, active, comm_solve, comm_rhs, grid_points, niter, dt; timeron = 0)

       npbversion="3.4.2"

       class = set_class(niter, grid_points)

       tsum = Array{Float64}(undef, t_last)
       t1 = Array{Float64}(undef, t_last)
       tming = Array{Float64}(undef, t_last)
       tmaxg = Array{Float64}(undef, t_last)

       t_recs = "total", "i/o", "rhs", "xsolve", "ysolve", "zsolve",
                "bpack", "exch", "xcomm", "ycomm", "zcomm",
                " totcomp", " totcomm"

       if (!active) @goto L999 end

#---------------------------------------------------------------------
#      Root node reads input file (if it exists) else takes
#      defaults from parameters
#---------------------------------------------------------------------
       if node == root

          @printf(stdout, "\n\n NAS Parallel Benchmarks 3.4 -- BT Benchmark\n\n", )

          @printf(stdout, " Size: %4ix%4ix%4i  (class %s)\n", grid_points[1], grid_points[2], grid_points[3], class)
          @printf(stdout, " Iterations: %4i    dt: %11.7F\n", niter, dt)
          @printf(stdout, " Total number of processes: %6i\n", total_nodes)
          if (no_nodes != total_nodes) 
             @printf(stdout, " WARNING: Number of processes is not a square number (%0i active)\n", no_nodes) 
          end
          println(stdout)
       end

       MAX_CELL_DIM, IMAX, JMAX, KMAX, BUF_SIZE, 
              cell_coord, cell_low, cell_high, cell_size, cell_start, cell_end, slice,
              forcing, u, rhs, lhsc, backsub_info, in_buffer, out_buffer, 
              #cv, rhon, rhos, rhoq, 
              cuf, q, ue, buf, fjac, njac, lhsa, lhsb, us, vs, ws, qs, rho_i, square = alloc_space(maxcells, grid_points[1])

       ncells = make_set(node, no_nodes, grid_points, cell_coord, cell_low, cell_high, cell_size, slice)

       for c = 1:maxcells
          if (cell_size[1, c] > IMAX) || (cell_size[2, c] > JMAX) ||(cell_size[3, c] > KMAX)
             println(stdout, node, c, view(cell_size, 1:3, c)...)
             println(stdout, " Problem size too big for compiled array sizes")
             @goto L999
          end
       end

       for i = 1:t_last
          timer_clear(i)
       end

       dnxm1, dnym1, dnzm1 , 
             tx1, tx2, tx3, 
             ty1, ty2, ty3, 
             tz1, tz2, tz3, 
             dttx1 , dttx2, 
             dtty1, dtty2, 
             dttz1, dttz2, 
             c2dttx1, c2dtty1, c2dttz1, dtdssp, 
             comz1, comz4, comz5, comz6, 
             c3c4tx3, c3c4ty3, c3c4tz3, 
             dx1tx1, dx2tx1, dx3tx1, dx4tx1, dx5tx1,
             dy1ty1, dy2ty1, dy3ty1, dy4ty1, dy5ty1, 
             dz1tz1, dz2tz1, dz3tz1, dz4tz1, dz5tz1 , 
             xxcon1, xxcon2, xxcon3, xxcon4, xxcon5,
             yycon1, yycon2, yycon3, yycon4, yycon5,
             zzcon1, zzcon2, zzcon3, zzcon4, zzcon5 = set_constants(dt, grid_points)

       initialize(IMAX, JMAX, KMAX, ncells, u, cell_low, cell_high, cell_size, slice, dnxm1, dnym1, dnzm1)

       setup_btio()

       lhsinit(ncells, lhsc, cell_coord, cell_start, cell_end, cell_size)

       exact_rhs(ncells, forcing, ue, buf, cuf, q, cell_low, cell_start, cell_end, cell_size, dssp, 
                   tx2, ty2, tz2, 
                   dnxm1, dnym1, dnzm1, 
                   xxcon1, xxcon2, xxcon3, xxcon4, xxcon5, 
                   yycon1, yycon2, yycon3, yycon4, yycon5, 
                   zzcon1, zzcon2, zzcon3, zzcon4, zzcon5, 
                   dx1tx1, dx2tx1, dx3tx1, dx4tx1, dx5tx1,
                   dy1ty1, dy2ty1, dy3ty1, dy4ty1, dy5ty1,
                   dz1tz1, dz2tz1, dz3tz1, dz4tz1, dz5tz1)

       ss, sr, b_size = compute_buffer_size(no_nodes, ncells, cell_coord, cell_size, 5)

       utmp = OffsetArray(zeros(Float64, 6, JMAX+4), 1:6, -2:JMAX+1)

#---------------------------------------------------------------------
#      do one time step to touch all code, and reinitialize
#---------------------------------------------------------------------
       adi(ss, sr, b_size,
            MAX_CELL_DIM,
            IMAX,
            JMAX,
            KMAX,
            cell_coord,
            cell_size,
            cell_start,
            cell_end,
            slice,
            forcing,           
            u,
            rhs,
            lhsc,
            backsub_info,
            in_buffer,
            out_buffer,
            fjac,
            njac,
            lhsa,
            lhsb,
            us,
            vs,
            ws,
            qs,
            rho_i,
            square,
            dt,
            timeron,
            Val(ncells),
            tx1,
            tx2,
            ty1,
            ty2,
            tz1,
            tz2,
            dx1tx1,
            dx2tx1,
            dx3tx1,
            dx4tx1,
            dx5tx1,
            dy1ty1,
            dy2ty1,
            dy3ty1,
            dy4ty1,
            dy5ty1,
            dz1tz1,
            dz2tz1,
            dz3tz1,
            dz4tz1,
            dz5tz1,
            xxcon2,
            xxcon3,
            xxcon4,
            xxcon5,
            yycon2,
            yycon3,
            yycon4,
            yycon5,
            zzcon2,
            zzcon3,
            zzcon4,
            zzcon5,
            Val(no_nodes), 
            comm_solve,
            comm_rhs,
            predecessor,
            successor,
            utmp,
       )
       initialize(IMAX, JMAX, KMAX, ncells, u, cell_low, cell_high, cell_size, slice, dnxm1, dnym1, dnzm1)

#---------------------------------------------------------------------
#      Synchronize before placing time stamp
#---------------------------------------------------------------------
       for i = 1:t_last
          timer_clear(i)
       end

       MPI.Barrier(comm_setup)

       timer_start(1)

       timer_clear(64)

       for STEP = 1:niter

          if node == root
             if mod(STEP, 20) == 0 || STEP == niter || STEP == 1
                @printf(stdout, " Time step %4i\n", STEP)
             end
          end

          adi(ss, sr, b_size,
               MAX_CELL_DIM,
               IMAX,
               JMAX,
               KMAX,
               cell_coord,
               cell_size,
               cell_start,
               cell_end,
               slice,
               forcing,           
               u,
               rhs,
               lhsc,
               backsub_info,
               in_buffer,
               out_buffer,
               fjac,
               njac,
               lhsa,
               lhsb,
               us,
               vs,
               ws,
               qs,
               rho_i,
               square,
               dt,
               timeron,
               Val(ncells),
               tx1,
               tx2,
               ty1,
               ty2,
               tz1,
               tz2,
               dx1tx1,
               dx2tx1,
               dx3tx1,
               dx4tx1,
               dx5tx1,
               dy1ty1,
               dy2ty1,
               dy3ty1,
               dy4ty1,
               dy5ty1,
               dz1tz1,
               dz2tz1,
               dz3tz1,
               dz4tz1,
               dz5tz1,
               xxcon2,
               xxcon3,
               xxcon4,
               xxcon5,
               yycon2,
               yycon3,
               yycon4,
               yycon5,
               zzcon2,
               zzcon3,
               zzcon4,
               zzcon5,
               Val(no_nodes), 
               comm_solve,
               comm_rhs,
               predecessor,
               successor,
               utmp,
  )

       end

       timer_start(2)
       btio_cleanup()
       timer_stop(2)

       timer_stop(1)
       t = timer_read(1)
       t1[1] = timer_read(t_enorm)

       timer_clear(t_enorm)
       verified = verify(no_nodes, node, 
                         ncells, 
                         class, 
                         grid_points,
                         cell_coord, cell_low, cell_high, cell_start, cell_end, cell_size,
                         forcing,        
                         u,
                         rhs,
                         in_buffer,
                         out_buffer,
                         us,
                         vs,
                         ws,
                         qs,
                         rho_i,
                         square,
                         dt,
                         ss,
                         sr,
                         b_size,
                         tx2, ty2, tz2,
                         dx1tx1, dx2tx1, dx3tx1, dx4tx1, dx5tx1,
                         dy1ty1, dy2ty1, dy3ty1, dy4ty1, dy5ty1,
                         dz1tz1, dz2tz1, dz3tz1, dz4tz1, dz5tz1,
                         dnxm1, dnym1, dnzm1,
                         xxcon2, xxcon3, xxcon4, xxcon5,
                         yycon2, yycon3, yycon4, yycon5,
                         zzcon2, zzcon3, zzcon4, zzcon5,
                         timeron,
                         comm_setup,
                         comm_rhs)

       tmax = MPI.Reduce(t, MPI.MAX, root, comm_setup)

       if node == root
          n3 = float(grid_points[1])*grid_points[2]*grid_points[3]
          navg = (grid_points[1]+grid_points[2]+grid_points[3])/3.0e0
          if tmax != 0.
             mflops = 1.0e-6*float(niter)*(3478.8*n3-17655.7*navg^2+28023.7*navg)/tmax
          else
             mflops = 0.0e0
          end

          print_results("BT", class, grid_points[1],
                              grid_points[2], grid_points[3], niter, no_nodes,
                              total_nodes, tmax, mflops, "          floating point",
                              verified, npbversion)
       end

       ttt = timer_read(64)
       @printf(stdout, " Time x_solve_cell =             %12.2F\n", ttt)

       if timeron
         for i = 1:t_zcomm
            t1[i] = timer_read(i)
         end
         t1[t_xsolve] = t1[t_xsolve] - t1[t_xcomm]
         t1[t_ysolve] = t1[t_ysolve] - t1[t_ycomm]
         t1[t_zsolve] = t1[t_zsolve] - t1[t_zcomm]
         t1[t_comm] = t1[t_xcomm]+t1[t_ycomm]+t1[t_zcomm]+t1[t_exch]
         t1[t_comp] = t1[t_total] - t1[t_comm]

         tsum = MPI.Reduce(t1, MPI.SUM, 0, comm_setup)
         tming = MPI.Reduce(t1, MPI.MIN, 0, comm_setup)
         tmaxg = MPI.Reduce(t1, MPI.MAX, 0, comm_setup)

         if node == 0
            @printf(stdout, " nprocs =%6i           minimum     maximum     average\n", no_nodes)
            for i = 1:t_last
               tsum[i] = tsum[i] / no_nodes
               @printf(stdout, " timer %2i(%8s) :  %10.4F  %10.4F  %10.4F\n", i, t_recs[i], tming[i], tmaxg[i], tsum[i])
            end
         end
       end

       @label L999
       MPI.Barrier(MPI.COMM_WORLD)
       #MPI.Finalize()

       return nothing
end

