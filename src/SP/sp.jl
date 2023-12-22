#-------------------------------------------------------------------------!
#                                                                         !
#        N  A  S     P A R A L L E L     B E N C H M A R K S  3.4         !
#                                                                         !
#                                   S P                                   !
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
# Authors: R. F. Van der Wijngaart
#          W. Saphir
#---------------------------------------------------------------------

tsum = Array{Float64}(undef, t_last+2)
t1 = Array{Float64}(undef, t_last+2)
tming = Array{Float64}(undef, t_last+2)
tmaxg = Array{Float64}(undef, t_last+2)

t_recs = ["total", "rhs", "xsolve", "ysolve", "zsolve",
          "bpack", "exch", "xcomm", "ycomm", "zcomm",
          " totcomp", " totcomm"]

#---------------------------------------------------------------------
function go()
#---------------------------------------------------------------------

       setup_mpi()

       if (!active) @goto L999 end

#---------------------------------------------------------------------
#      Root node reads input file (if it exists) else takes
#      defaults from parameters
#---------------------------------------------------------------------
       if node == root

          @printf(stdout, "\n\n NAS Parallel Benchmarks 3.4 -- SP Benchmark\n\n", )

          global timeron = check_timer_flag()

          fstatus = isfile("inputsp.data") ? 0 : 1
#
          if fstatus == 0
            @printf(stdout, " Reading from input file inputsp.data\n", )
            f = open("inputsp.data","r")
            global niter = parse(Int, readline(f))
            global dt = parse(Int, readline(f))
            grid_points[1] = parse(Int, readline(f))
            grid_points[2] = parse(Int, readline(f))
            grid_points[3] = parse(Int, readline(f))
            close(f)
          else
            @printf(stdout, " No input file inputsp.data. Using compiled defaults\n", )
            global niter = niter_default
            global dt    = dt_default
            grid_points[1] = problem_size
            grid_points[2] = problem_size
            grid_points[3] = problem_size
          end

          class = set_class(niter)

          @printf(stdout, " Size: %4ix%4ix%4i  (class %s)\n", grid_points[1], grid_points[2], grid_points[3],  class)
          @printf(stdout, " Iterations: %4i    dt: %11.7F\n", niter, dt)
          @printf(stdout, " Total number of processes: %6i\n", total_nodes)
          if (no_nodes != total_nodes) 
            @printf(stdout, " WARNING: Number of processes is not a square number (%0i active)\n", no_nodes) 
          end
          println(stdout)

       else
         global niter = -1
         global dt = -1
         global timeron = -1
         class = "U"
       end

       niter = MPI.bcast(niter, comm_setup; root=root)

       dt = MPI.bcast(dt, comm_setup; root=root)

       grid_points_0 = MPI.bcast(grid_points, comm_setup; root=root)
       grid_points[1] = grid_points_0[1]
       grid_points[2] = grid_points_0[2]
       grid_points[3] = grid_points_0[3]

       timeron = MPI.bcast(timeron, comm_setup; root=root)

       alloc_space()

       make_set()

       for c = 1:ncells
          if (cell_size[1, c] > IMAX) ||(cell_size[2, c] > JMAX) ||(cell_size[3, c] > KMAX)
             println(stdout, node, c, view(cell_size, 1:3, c)...)
             println(stdout, " Problem size too big for compiled array sizes")
             @goto L999
          end
       end

       for i = 1:t_last
          timer_clear(i)
       end

       set_constants()

       initialize()

       lhsinit()

       exact_rhs()

       compute_buffer_size(5)

#---------------------------------------------------------------------
#      do one time step to touch all code, and reinitialize
#---------------------------------------------------------------------
       adi()
       initialize()

#---------------------------------------------------------------------
#      Synchronize before placing time stamp
#---------------------------------------------------------------------
       for i = 1:t_last
          timer_clear(i)
       end

       MPI.Barrier(comm_setup)

       timer_clear(1)
       timer_start(1)

       for STEP = 1:niter

          if node == root
             if mod(STEP, 20) == 0 || STEP == 1
                @printf(stdout, " Time step %4i\n", STEP)
              end
          end

          adi()

       end

       timer_stop(1)
       t = timer_read(1)

       verified = verify(class)

       tmax = MPI.Reduce(t, MPI.MAX, root, comm_setup)

       if node == root
          if tmax != 0.
             n3 = float(grid_points[1])*grid_points[2]*grid_points[3]
             t = (grid_points[1]+grid_points[2]+grid_points[3])/3.0e0
             mflops = 1.0e-6*float(niter)*(881.174*n3-4683.91* t^2+11484.5* t-19272.4) / tmax
          else
             mflops = 0.0e0
          end

         print_results("SP", class, grid_points[1],
           grid_points[2], grid_points[3], niter, no_nodes,
           total_nodes, tmax, mflops, "          floating point",
           verified, npbversion)
       end

       if (!timeron) @goto L999 end

       for i = 1:t_last
          t1[i] = timer_read(i)
       end
       t1[t_xsolve] = t1[t_xsolve] - t1[t_xcomm]
       t1[t_ysolve] = t1[t_ysolve] - t1[t_ycomm]
       t1[t_zsolve] = t1[t_zsolve] - t1[t_zcomm]
       t1[t_last+2] = t1[t_xcomm]+t1[t_ycomm]+t1[t_zcomm]+t1[t_exch]
       t1[t_last+1] = t1[t_total]  - t1[t_last+2]

       tsum = MPI.Reduce(t1, MPI.SUM, 0, comm_setup)
       tming = MPI.Reduce(t1, MPI.MIN, 0, comm_setup)
       tmaxg = MPI.Reduce(t1, MPI.MAX, 0, comm_setup)

       if node == 0
          @printf(stdout, " nprocs =%6i           minimum     maximum     average\n", no_nodes)
          for i = 1:t_last+2
             tsum[i] = tsum[i] / no_nodes
             @printf(stdout, " timer %2i(%8s) :  %10.4F  %10.4F  %10.4F\n", i, t_recs[i], tming[i], tmaxg[i], tsum[i])
          end
       end

       @label L999
       MPI.Barrier(MPI.COMM_WORLD)
       MPI.Finalize()

       return nothing
end
