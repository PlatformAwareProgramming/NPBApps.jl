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


function go(class::CLASS)

   setup_mpi()

   problem_size = bt_class[class].problem_size
   
   niter = bt_class[class].niter
   dt    = bt_class[class].dt

   grid_points = zeros(Integer, 3)
   grid_points[1] = problem_size
   grid_points[2] = problem_size
   grid_points[3] = problem_size

   go(grid_points, niter, dt)

end

function go(params_file::String)

   setup_mpi()

   if node == root

      fstatus = isfile(params_file) ? 0 : 1
   #
      grid_points = zeros(Integer, 3)
      if fstatus == 0
         @printf(stdout, " Reading from input file params_file\n", )
         f = open(params_file,"r")
         niter = parse(Int, readline(f))
         dt = parse(Int, readline(f))
         grid_points[1] = parse(Int, readline(f))
         grid_points[2] = parse(Int, readline(f))
         grid_points[3] = parse(Int, readline(f))
         close(f)
      else
         @printf(stdout, " No input file params_file. Using defaults (class S) \n", )
         problem_size =  class[S].problem_size
         niter = class[S].niter
         dt    = class[S].dt
         grid_points[1] = problem_size
         grid_points[2] = problem_size
         grid_points[3] = problem_size
      end
   else
      niter = -1
      dt = -1
      class = CLASS_UNDEFINED
   end

   niter = MPI.bcast(niter, comm_setup; root=root)
   dt = MPI.bcast(dt, comm_setup; root=root)

   grid_points_0 = MPI.bcast(grid_points, comm_setup; root=root)
   grid_points[1] = grid_points_0[1]
   grid_points[2] = grid_points_0[2]
   grid_points[3] = grid_points_0[3]

   perform(grid_points, niter, dt)
end

function go()
   go("inputbt.data")
end

function go(grid_points, niter, dt)

   setup_mpi()

   perform(grid_points, niter, dt)

end

function perform(grid_points, niter, dt)

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

          global timeron = check_timer_flag()

          @printf(stdout, " Size: %4ix%4ix%4i  (class %s)\n", grid_points[1], grid_points[2], grid_points[3], class)
          @printf(stdout, " Iterations: %4i    dt: %11.7F\n", niter, dt)
          @printf(stdout, " Total number of processes: %6i\n", total_nodes)
          if (no_nodes != total_nodes) 
             @printf(stdout, " WARNING: Number of processes is not a square number (%0i active)\n", no_nodes) 
          end
          println(stdout)
       else
            global timeron = -1
       end

       timeron = MPI.bcast(timeron, comm_setup; root=root)

       alloc_space(grid_points[1])

       make_set(grid_points)

       for c = 1:maxcells
          if (cell_size[z][1, c] > IMAX) || (cell_size[z][2, c] > JMAX) ||(cell_size[z][3, c] > KMAX)
             println(stdout, node, c, view(cell_size, 1:3, c)...)
             println(stdout, " Problem size too big for compiled array sizes")
             @goto L999
          end
       end

       for i = 1:t_last
          timer_clear(i)
       end

       set_constants(dt, grid_points)

       initialize()

       setup_btio()
       idump = 0

       lhsinit()

       exact_rhs()

       compute_buffer_size(5)

       if (no_nodes > 1)
         ss = SA[start_send_east::Int start_send_west::Int start_send_north::Int start_send_south::Int start_send_top::Int start_send_bottom::Int]
         sr = SA[start_recv_east::Int start_recv_west::Int start_recv_north::Int start_recv_south::Int start_recv_top::Int start_recv_bottom::Int]
         b_size = SA[east_size::Int west_size::Int north_size::Int south_size::Int top_size::Int bottom_size::Int]
       else
         ss = nothing
         sr = nothing
         b_size = nothing
       end

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
       initialize()

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
       verified = verify(class, grid_points, dt, ss, sr, b_size)

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

       if (!timeron) @goto L999 end

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

       @label L999
       MPI.Barrier(MPI.COMM_WORLD)
       MPI.Finalize()

       return nothing
end

