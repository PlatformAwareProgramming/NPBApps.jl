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
# Authors: R. F. Van der Wijngaart
#          T. Harris
#          M. Yarrow
#---------------------------------------------------------------------

function perform(clusterid_, clusters,  niter,  dt,  ratio,  x_zones,  y_zones,  gx_size,  gy_size,  gz_size,  problem_size,  nxmax,  nx,  ny,  nz,  proc_num_zones,  proc_zone_id_, zone_proc_id_, iz_west,  iz_east,  iz_north,  iz_south, itimer_=false, npb_verbose_=0)
       global clusterid = clusterid_
       global num_clusters = length(clusters) 
       global itimer = itimer_
       global timeron = itimer > 0
       global npb_verbose = npb_verbose_
       global max_zones = x_zones * y_zones
       global proc_zone_id = proc_zone_id_
       global zone_proc_id = zone_proc_id_
       num_zones = max_zones
       
       setup_mpi(proc_num_zones) 

       #GC.enable(false)

       @info "$clusterid/$node: NUM_THREADS = $(Threads.nthreads()) ---- number of zones = $proc_num_zones --- $(MPI.Query_thread()) "

       if (!active) @goto L999 end

       class = set_class(niter, x_zones, y_zones, gx_size, gy_size, gz_size)

       alloc_field_space_zones(proc_num_zones)
       
       total_size = 0
       for iz = 1:proc_num_zones

          zone = proc_zone_id[iz]

          total_size += alloc_field_space(iz, [nx[zone], ny[zone], nz[zone]], problem_size)
          make_set(iz, [nx[zone], ny[zone], nz[zone]])

          for c = 1:ncells
             if (cell_size[iz][1, c] > IMAX[iz]) ||(cell_size[iz][2, c] > JMAX[iz]) ||(cell_size[iz][3, c] > KMAX[iz])
                println(stdout, node, c, view(cell_size[iz], 1:3, c)...)
                println(stdout, " Problem size too big for compiled array sizes")
                #@goto L999
             end
         end
      end

      @info "$clusterid/$node: TOTAL SIZE = $(total_size) bytes"

      for i = 1:t_last
          timer_clear(i)
       end

       ss = Array{SVector{6,Int}}(undef, proc_num_zones)
       sr = Array{SVector{6,Int}}(undef, proc_num_zones)
       b_size = Array{SVector{6,Int}}(undef, proc_num_zones)

       compute_buffer_size_initial(proc_num_zones)
       set_constants(dt, gx_size, gy_size, gz_size, x_zones, y_zones) 

       Threads.@threads for iz = 1:proc_num_zones         
         initialize(iz) 
         lhsinit(iz) 
         exact_rhs(iz) 
         compute_buffer_size(iz, 5)

         ss[iz] = SA[start_send_east[iz]::Int start_send_west[iz]::Int start_send_north[iz]::Int start_send_south[iz]::Int start_send_top[iz]::Int start_send_bottom[iz]::Int]
         sr[iz] = SA[start_recv_east[iz]::Int start_recv_west[iz]::Int start_recv_north[iz]::Int start_recv_south[iz]::Int start_recv_top[iz]::Int start_recv_bottom[iz]::Int]
         b_size[iz] = SA[east_size[iz]::Int west_size[iz]::Int north_size[iz]::Int south_size[iz]::Int top_size[iz]::Int bottom_size[iz]::Int]
      end

      requests = Array{Array{MPI.Request}}(undef,proc_num_zones)
       s = Array{Array{FloatType}}(undef,proc_num_zones)
       utmp = Array{OffsetArray{FloatType, 2, Array{FloatType, 2}}}(undef,proc_num_zones)
       send_id = Array{Ref{MPI.Request}}(undef,proc_num_zones)
       recv_id = Array{Ref{MPI.Request}}(undef,proc_num_zones)
       for iz = 1:proc_num_zones
         requests[iz] = Array{MPI.Request}(undef,12)
         s[iz] = Array{FloatType}(undef,5)
         utmp[iz] = OffsetArray(zeros(FloatType, 6, JMAX[iz]+4), 1:6, -2:JMAX[iz]+1)
         send_id[iz] = Ref{MPI.Request}(MPI.REQUEST_NULL)
         recv_id[iz] = Ref{MPI.Request}(MPI.REQUEST_NULL)
      end

#---------------------------------------------------------------------
#      do one time step to touch all code, and reinitialize
#---------------------------------------------------------------------

      if (no_nodes > 1 && num_clusters > 1) || proc_num_zones > 1
            exch_qbc(Val(ncells),
                     proc_num_zones,
                     zone_proc_id,
                     proc_zone_id,
                     iz_west,
                     iz_east,
                     iz_north,
                     iz_south,
                     cell_coord,
                     cell_size, 
                     cell_start,
                     cell_end,
                     cell_low,
                     cell_high,
                     successor,
                     predecessor,
                     u,
                     in_buffer,
                     out_buffer,
                     requests,
                     ss,
                     sr,
                     b_size,
                     comm_exch,
                     timeron,)  
       end

       Threads.@threads for iz = 1:proc_num_zones
            adi(iz, ss[iz], 
                  sr[iz], 
                  b_size[iz],
                  MAX_CELL_DIM[iz],
                  IMAX[iz],
                  JMAX[iz],
                  KMAX[iz],
                  cell_coord[iz],
                  cell_size[iz],
                  cell_start[iz],
                  cell_end[iz],
                  slice[iz],
                  forcing[iz],           
                  u[iz],
                  rhs[iz],
                  lhsc[iz],
                  backsub_info[iz],
                  in_buffer[iz],
                  out_buffer[iz],
                  fjac[iz],
                  njac[iz],
                  lhsa[iz],
                  lhsb[iz],
                  us[iz],
                  vs[iz],
                  ws[iz],
                  qs[iz],
                  rho_i[iz],
                  square[iz],
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
                  comm_solve[iz],
                  comm_rhs[iz],
                  predecessor[iz],
                  successor[iz],
                  utmp[iz],
                  requests[iz],
                  send_id[iz],
                  recv_id[iz],
               )
       end


       #@goto L999

       for iz = 1:proc_num_zones
           initialize(iz)
       end


#---------------------------------------------------------------------
#      Synchronize before placing time stamp
#---------------------------------------------------------------------
       for i = 1:t_last
          timer_clear(i)
       end

       MPI.Barrier(comm_setup)

       timer_clear(1)
       timer_start(1)

       @label L997

       Q = 1

       timer_clear(64); t_64 = 0.0; t_64s = 0.0
       timer_clear(63); t_63 = 0.0; t_63s = 0.0

       for STEP = 1:niter
          #GC.gc()
          if node == root
            Q = STEP > 1 && t_64 < 5.0 ? ceil(5.0 / t_64) : Q

            if mod(STEP, Q) == 0 || STEP == 1
               @printf(stdout, "%2i: Time step %4i  -- %12.2F  -- %12.2F --- %4i \n", clusterid, STEP, t_63s/(STEP-1), t_64s/(STEP-1), Q)
            end
          end

          timer_start(63)
 
          if (no_nodes > 1 && num_clusters > 1) || proc_num_zones > 1
            exch_qbc(Val(ncells),
                     proc_num_zones,
                     zone_proc_id,
                     proc_zone_id,
                     iz_west,
                     iz_east,
                     iz_north,
                     iz_south,
                     cell_coord,
                     cell_size,  
                     cell_start,
                     cell_end,
                     cell_low,
                     cell_high,
                     successor,
                     predecessor,
                     u,
                     in_buffer,
                     out_buffer,
                     requests,
                     ss,
                     sr,
                     b_size,
                     comm_exch,
                     timeron,) 
          end
 
          t_63 = timer_stop(63); t_63s += t_63

          timer_start(64)

          Threads.@threads for iz = 1:proc_num_zones                          
               adi(iz, ss[iz], 
                  sr[iz], 
                  b_size[iz],
                  MAX_CELL_DIM[iz],
                  IMAX[iz],
                  JMAX[iz],
                  KMAX[iz],
                  cell_coord[iz],
                  cell_size[iz],
                  cell_start[iz], 
                  cell_end[iz], 
                   slice[iz],
                  forcing[iz],           
                  u[iz],
                  rhs[iz],
                  lhsc[iz],
                  backsub_info[iz],
                  in_buffer[iz],
                  out_buffer[iz],
                  fjac[iz],
                  njac[iz],
                  lhsa[iz],
                  lhsb[iz],
                  us[iz],
                  vs[iz],
                  ws[iz],
                  qs[iz],
                  rho_i[iz],
                  square[iz],
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
                  comm_solve[iz],
                  comm_rhs[iz],
                  predecessor[iz], 
                  successor[iz],
                  utmp[iz],
                  requests[iz],
                  send_id[iz],
                  recv_id[iz],
               )
           end

           t_64 = timer_stop(64); t_64s += t_64

       end

       timer_stop(1)
       t = timer_read(1)

#---------------------------------------------------------------------
#      perform verification and print results
#---------------------------------------------------------------------

       verify(dt, ss, sr, b_size, proc_num_zones, proc_zone_id, rho_i, us, vs, ws, qs, square, rhs, forcing, u, nx, ny, nz, requests)

       tmax = MPI.Reduce(t, MPI.MAX, root, comm_setup)
       if node == root
         remotecall(reportMaxTimeNode, 1, tmax; role=:worker)
       end

#---------------------------------------------------------------------
#      More timers
#---------------------------------------------------------------------
       if (!timeron) @goto L999 end

       for i = 1:t_last
          trecs[i] = timer_read(i)
       end

       tsum = MPI.Reduce(trecs, MPI.SUM, root, comm_setup) 
       tming = MPI.Reduce(trecs, MPI.MIN, root, comm_setup)
       tmaxg = MPI.Reduce(trecs, MPI.MAX, root, comm_setup)

       if node == root
         tavg = tsum ./ no_nodes
         remotecall(reportTimersNode, 1, tavg, tming, tmaxg; role=:worker)
       end

       @label L999
       MPI.Barrier(MPI.COMM_WORLD)
       MPI.Finalize()

       return nothing
end











