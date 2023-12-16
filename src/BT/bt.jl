using FortranFiles
using OffsetArrays
using Parameters
using Printf

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

#---------------------------------------------------------------------
       function MPBT()
#---------------------------------------------------------------------

#       use bt_data
#       use mpinpb

#       implicit none

#       integer i, niter, STEP, c, ERROR, fstatus
#       DOUBLEPRECISION navg, mflops, mbytes, n3

#       external timer_read
#       DOUBLEPRECISION t, tmax, iorate[2], tpc, timer_read
#       logical verified
#       character class, cbuff*40
#       DOUBLEPRECISION t1[t_last], tsum[t_last],  
#                        tming[t_last], tmaxg[t_last]
#       character        t_recs[t_last]

#       integer wr_interval

            t_recs = "total", "i/o", "rhs", "xsolve", "ysolve", "zsolve",
                   "bpack", "exch", "xcomm", "ycomm", "zcomm",
                   " totcomp", " totcomm"

       setup_mpi
       if (!active) @goto L999 end

#---------------------------------------------------------------------
#      Root node reads input file (if it exists) else takes
#      defaults from parameters
#---------------------------------------------------------------------
       if node == root

          @printf(stdout, "\n\n NAS Parallel Benchmarks 3.4 -- BT Benchmark\n\n", )

          check_timer_flag( timeron )

          OPEN(unit = 2, file = "inputbt.data", status = "old", iostat = fstatus)
#
          rd_interval = 0
          if fstatus == 0
            @printf(stdout, " Reading from input file inputbt.data\n", )
# 233        format(' Reading from input file inputbt.data')
            READ(2, niter)
            READ(2, dt)
            READ(2, predecessor[1], predecessor[2], predecessor[3])
            if iotype != 0
                READ(2, "%s", cbuff)
                READ(cbuff, i, wr_interval, rd_interval)
                if (i != 0) rd_interval = 0 end
                if (wr_interval <= 0) wr_interval = wr_default end
            end
            if iotype == 1
                READ(2, collbuf_nodes, collbuf_size)
                println(stdout, "collbuf_nodes ", collbuf_nodes)
                println(stdout, "collbuf_size  ", collbuf_size)
            end
            CLOSE(2)
          else
            @printf(stdout, " No input file inputbt.data. Using compiled defaults\n", )
            niter = niter_default
            dt    = dt_default
            predecessor[1] = problem_size
            predecessor[2] = problem_size
            predecessor[3] = problem_size
            wr_interval = wr_default
            if iotype == 1
#             set number of nodes involved in collective buffering to 4,
#             unless total number of nodes is smaller than that.
#             set buffer size for collective buffering to 1MB per node
#             collbuf_nodes = min(4,no_nodes)
#             set default to No-File-Hints with a value of 0
              collbuf_nodes = 0
              collbuf_size = 1000000
            end
          end
# 234      format(' No input file inputbt.data. Using compiled defaults')

          set_class(niter, class)

          @printf(stdout, " Size: %4ix%4ix%4i  (class %s)\n", predecessor[1], predecessor[2], predecessor[3],
                         class)
          @printf(stdout, " Iterations: %4i    dt: %11.7F\n", niter, dt)
          @printf(stdout, " Total number of processes: %6i\n", total_nodes)
          if (no_nodes != total_nodes) @printf(stdout, " WARNING: Number of processes is not a square number (%0i active)\n", no_nodes) end
          println(stdout, )

          if (iotype == 1) @printf(stdout, " BTIO -- %s write interval: %3i\n\n", "FULL MPI-IO", wr_interval) end
          if (iotype == 2) @printf(stdout, " BTIO -- %s write interval: %3i\n\n", "SIMPLE MPI-IO", wr_interval) end
          if (iotype == 3) @printf(stdout, " BTIO -- %s write interval: %3i\n\n", "EPIO", wr_interval) end
          if (iotype == 4) @printf(stdout, " BTIO -- %s write interval: %3i\n\n", "FORTRAN IO", wr_interval) end

# 1000 format(//, ' NAS Parallel Benchmarks 3.4 -- BT Benchmark',/)
# 1001     format(' Size: ', i4, 'x', i4, 'x', i4, '  (class ', a, ')' )
# 1002     format(' Iterations: ', i4, '    dt: ', F11.7)
# 1003     format(' Total number of processes: ', i6)
# 1004     format(' WARNING: Number of processes is not a square number',  		                 ' (', i0, ' active)')
# 1006     format(' BTIO -- ', A, ' write interval: ', i3 /)

       end

       mpi_bcast(niter, 1, MPI_INTEGER,
                      root, comm_setup, ERROR)

       mpi_bcast(dt, 1, dp_type,
                      root, comm_setup, ERROR)

       mpi_bcast(predecessor[1], 3, MPI_INTEGER,
                      root, comm_setup, ERROR)

       mpi_bcast(wr_interval, 1, MPI_INTEGER,
                      root, comm_setup, ERROR)

       mpi_bcast(rd_interval, 1, MPI_INTEGER,
                      root, comm_setup, ERROR)

       mpi_bcast(timeron, 1, MPI_LOGICAL,
                      root, comm_setup, ERROR)

       alloc_space

       make_set

       for c = 1:maxcells
          if (cell_size[1, c] > IMAX) ||(
               cell_size[2, c] > JMAX) ||(
               cell_size[3, c] > KMAX)
             println(stdout, node, c, view(cell_size, 1:3, c)...)
             println(stdout, " Problem size too big for compiled array sizes")
              @goto L999
          end
       end

       for i = 1:t_last
          timer_clear(i)
       end

       set_constants

       initialize

       setup_btio
       idump = 0

       lhsinit

       exact_rhs

       compute_buffer_size(5)

#---------------------------------------------------------------------
#      do one time step to touch all code, and reinitialize
#---------------------------------------------------------------------
       adi
       initialize

#---------------------------------------------------------------------
#      Synchronize before placing time stamp
#---------------------------------------------------------------------
       for i = 1:t_last
          timer_clear(i)
       end
       mpi_barrier(comm_setup, ERROR)

       timer_start(1)

       for STEP = 1:niter

          if node == root
             if mod(STEP, 20) == 0 || STEP == niter ||
                 STEP == 1
                @printf(stdout, " Time step %4i\n", STEP)
# 200            format(' Time step ', i4)
             end
          end

          adi

          if iotype != 0
              if mod(STEP, wr_interval) == 0 || STEP == niter
                  if node == root
                      println(stdout, "Writing data set, time step", STEP)
                  end
                  if STEP == niter && rd_interval > 1
                      rd_interval = 1
                  end
                  timer_start(2)
                  output_timestep
                  timer_stop(2)
                  idump = idump + 1
              end
          end
       end

       timer_start(2)
       btio_cleanup
       timer_stop(2)

       timer_stop(1)
       t = timer_read(1)
       t1[1] = timer_read(t_enorm)

       timer_clear(t_enorm)
       verify(class, verified)

       mpi_reduce(t, tmax, 1,
                       dp_type, MPI_MAX,
                       root, comm_setup, ERROR)

       if iotype != 0
          n3 = 0.0e0
          for c = 1:ncells
             n3 = n3 + float(cell_size[1, c]) * cell_size[2, c] * cell_size[3, c]
          end
          mbytes = n3 * 40.0 * idump * 1.0e-6
          t1[2] = timer_read(t_enorm)
          for i = 1:2
             if i == 1
                t = timer_read(t_io)
             else
                t = timer_read(t_iov)
             end
             t = t - t1[i]                      # remove enorm time
             if (t != 0.0e0) t = mbytes / t end  # rate MB/s
             t1[i] = t
          end
          if (rd_interval > 0) t1[1] = t1[1] * 2 end
          mpi_reduce(t1, iorate, 2,
                          dp_type, MPI_SUM,
                          root, comm_setup, ERROR)
       end

       if node == root
          n3 = float(predecessor[1])*predecessor[2]*predecessor[3]
          navg = (predecessor[1]+predecessor[2]+predecessor[3])/3.0e0
          if tmax != 0.
             mflops = 1.0e-6*float(niter)*(
                      3478.8*n3-17655.7*navg^2+28023.7*navg)/
                       tmax
          else
             mflops = 0.0e0
          end

          if iotype != 0
             mbytes = n3 * 40.0 * idump * 1.0e-6
             for i = 1:2
                t1[i] = 0.0
                if (iorate[i] != 0.0e0) t1[i] = mbytes / iorate[i] end
             end
             if (rd_interval > 0) t1[1] = t1[1] * 2 end
             tpc = 0.0
             if (tmax != 0.) tpc = t1[1] * 100.0 / tmax end
             @printf(stdout, "\n BTIO -- statistics:\n   I/O timing in seconds   : %14.2F\n   I/O timing percentage   : %14.2F\n   I/O timing in verify    : %14.2F\n   Total data written (MB) : %14.2F\n   I/O data rate  (MB/sec) : %14.2F\n", t1[1], tpc, t1[2], mbytes, iorate[1])
# 1100        format(/' BTIO -- statistics:'/  		                     '   I/O timing in seconds   : ', f14.2/  		                     '   I/O timing percentage   : ', f14.2/  		                     '   I/O timing in verify    : ', f14.2/  		                     '   Total data written (MB) : ', f14.2/  		                     '   I/O data rate  (MB/sec) : ', f14.2)
          end

         print_results("BT", class, predecessor[1],
           predecessor[2], predecessor[3], niter, no_nodes,
           total_nodes, tmax, mflops, "          floating point",
           verified, npbversion, compiletime, cs1, cs2, cs3, cs4, cs5,
           cs6, "(none)")
       end

       if (!timeron) @goto L999 end

       for i = 1:t_zcomm
          t1[i] = timer_read(i)
       end
       t1[t_xsolve] = t1[t_xsolve] - t1[t_xcomm]
       t1[t_ysolve] = t1[t_ysolve] - t1[t_ycomm]
       t1[t_zsolve] = t1[t_zsolve] - t1[t_zcomm]
       t1[t_comm] = t1[t_xcomm]+t1[t_ycomm]+t1[t_zcomm]+t1[t_exch]
       t1[t_comp] = t1[t_total] - t1[t_comm]

       MPI_Reduce(t1, tsum,  t_last, dp_type, MPI_SUM,
                       0, comm_setup, ERROR)
       MPI_Reduce(t1, tming, t_last, dp_type, MPI_MIN,
                       0, comm_setup, ERROR)
       MPI_Reduce(t1, tmaxg, t_last, dp_type, MPI_MAX,
                       0, comm_setup, ERROR)

       if node == 0
          @printf(stdout, " nprocs =%6i           minimum     maximum     average\n", no_nodes)
          for i = 1:t_last
             tsum[i] = tsum[i] / no_nodes
             @printf(stdout, " timer %2i(%8s) :  %10.4F  %10.4F  %10.4F\n", i, t_recs[i], tming[i], tmaxg[i], tsum[i])
          end
       end
# 800   format(' nprocs =', i6, 11x, 'minimum', 5x, 'maximum',  		              5x, 'average')
# 810   format(' timer ', i2, '(', A8, ') :', 3(2x,f10.4))

       @label L999
       mpi_barrier(MPI_COMM_WORLD, ERROR)
       mpi_finalize(ERROR)

       return nothing
       end

