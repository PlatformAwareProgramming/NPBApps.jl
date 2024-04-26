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






function perform(clusterid_, clusters, niter, dt, ratio, x_zones, y_zones, gx_size, gy_size, gz_size, problem_size, nxmax, nx, ny, nz, proc_num_zones, proc_zone_id, zone_proc_id, iz_west, iz_east, iz_north, iz_south, itimer_=false, npb_verbose_=0)

       global clusterid = clusterid_
       global num_clusters = length(clusters) 
       global itimer = itimer_
       global timeron = itimer > 0
       global npb_verbose = npb_verbose_
       global max_zones = x_zones * y_zones
       num_zones = max_zones
       
       @info "$clusterid/??: STEP 1 - $proc_num_zones" 

       setup_mpi(proc_num_zones) 

       @info "$clusterid/$node: STEP 2" 

       if (!active) @goto L999 end

       @info "$clusterid/$node: STEP 3" 

       class = set_class(niter, x_zones, y_zones, gx_size, gy_size, gz_size)

       @info "$clusterid/$node: STEP 4" 

       alloc_field_space_zones(proc_num_zones)
       
       @info "$clusterid/$node: STEP 5" 

       total_size = 0
       for iz = 1:proc_num_zones

          zone = proc_zone_id[iz]

          total_size += alloc_field_space(iz, [nx[zone], ny[zone], nz[zone]], problem_size) 
          make_set(iz, [nx[zone], ny[zone], nz[zone]])

          for c = 1:ncells
             if (cell_size[iz][1, c] > IMAX[iz]) ||(cell_size[iz][2, c] > JMAX[iz]) ||(cell_size[iz][3, c] > KMAX[iz])
                println(stdout, node, c, view(cell_size[iz], 1:3, c)...)
                println(stdout, " Problem size too big for compiled array sizes")
                @goto L999
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

       for iz = 1:proc_num_zones         
         set_constants(iz, dt, gx_size, gy_size, gz_size, x_zones, y_zones) 
         initialize(iz) 
         lhsinit(iz) 
         exact_rhs(iz) 
         compute_buffer_size(iz, 5)

         #if (no_nodes > 1)
            ss[iz] = SA[start_send_east[iz]::Int start_send_west[iz]::Int start_send_north[iz]::Int start_send_south[iz]::Int start_send_top[iz]::Int start_send_bottom[iz]::Int]
            sr[iz] = SA[start_recv_east[iz]::Int start_recv_west[iz]::Int start_recv_north[iz]::Int start_recv_south[iz]::Int start_recv_top[iz]::Int start_recv_bottom[iz]::Int]
            b_size[iz] = SA[east_size[iz]::Int west_size[iz]::Int north_size[iz]::Int south_size[iz]::Int top_size[iz]::Int bottom_size[iz]::Int]
         #else
         #   ss[iz] = SA[0 0 0 0 0 0]
         #   sr[iz] = SA[0 0 0 0 0 0]
         #   b_size[iz] = SA[0 0 0 0 0 0]
         #end
      end

       requests = Array{Array{MPI.Request}}(undef,proc_num_zones)
       s = Array{Array{Float32}}(undef,proc_num_zones)
       for iz = 1:proc_num_zones
         requests[iz] = Array{MPI.Request}(undef,12)
         for i = 1:12
             requests[iz][i] = MPI.REQUEST_NULL
         end
         s[iz] = Array{Float32}(undef,5)
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

       #@info "FINISHED !"
       #@goto L999

       Threads.@threads for zone = 1:proc_num_zones
         adi(zone, 
               Val(no_nodes),
               Val(ncells),
               slice[zone],
               cell_size[zone],
               cell_start[zone],
               cell_end[zone],
               cell_coord[zone],
               cell_low[zone],
               cell_high[zone],
               u[zone],
               rhs[zone],
               lhs[zone],
               rho_i[zone],
               us[zone],
               cv[zone],
               rhoq[zone],
               rhon[zone],
               rhos[zone],
               vs[zone],
               ws[zone],
               square[zone],
               qs[zone],
               ainv[zone],
               speed[zone],
               forcing[zone],
               dt,
               tx2,
               ty2,
               tz2,
               c1,
               c2,
               c1c2,
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
               dssp,
               con43,
               dx2, dx5, c3c4, c1c5, c2dtty1, dxmax, dx1, dttx1, dttx2,
               comz5, comz4, comz1, comz6,
               dy3, dy5, dy1, dymax, dtty1, dtty2, 
               dz4, dz5, dz1, dttz2, dttz1, dzmax,
               bt,
               successor[zone],
               predecessor[zone],
               in_buffer[zone],
               out_buffer[zone],
               comm_solve[zone],
               comm_rhs[zone],
               ss[zone],
               sr[zone],
               b_size[zone],
               s[zone],
               requests[zone]
            )
       end

       for zone = 1:proc_num_zones
           initialize(zone)
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

       @info "$clusterid/$node: NUM_THREADS = $(Threads.nthreads())"

       for STEP = 1:niter
          if node == root
            Q = STEP > 1 && t_64 < 5.0 ? ceil(5.0 / t_64) : Q

            if mod(STEP, Q) == 0 || STEP == 1
               @printf(stdout, "%2i: Time step %4i  -- %12.2F  -- %12.2F --- %4i \n", clusterid, STEP, t_63s/(STEP-1), t_64s/(STEP-1), Q)
            end
          end

          timer_start(64)
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

          Threads.@threads for zone = 1:proc_num_zones                          
              adi(zone, 
                  Val(no_nodes),
                  Val(ncells),
                  slice[zone],
                  cell_size[zone],
                  cell_start[zone],
                  cell_end[zone],
                  cell_coord[zone],
                  cell_low[zone],
                  cell_high[zone],
                  u[zone],
                  rhs[zone],
                  lhs[zone],
                  rho_i[zone],
                  us[zone],
                  cv[zone],
                  rhoq[zone],
                  rhon[zone],
                  rhos[zone],
                  vs[zone],
                  ws[zone],
                  square[zone],
                  qs[zone],
                  ainv[zone],
                  speed[zone],
                  forcing[zone],
                  dt,
                  tx2,
                  ty2,
                  tz2,
                  c1,
                  c2,
                  c1c2,
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
                  dssp,
                  con43,
                  dx2, dx5, c3c4, c1c5, c2dtty1, dxmax, dx1, dttx1, dttx2,
                  comz5, comz4, comz1, comz6,
                  dy3, dy5, dy1, dymax, dtty1, dtty2, 
                  dz4, dz5, dz1, dttz2, dttz1, dzmax,
                  bt,
                  successor[zone],
                  predecessor[zone],
                  in_buffer[zone],
                  out_buffer[zone],
                  comm_solve[zone],
                  comm_rhs[zone],
                  ss[zone],
                  sr[zone],
                  b_size[zone],
                  s[zone],
                  requests[zone]
            )
         end

         t_64 = timer_stop(64); t_64s += t_64

       end

       timer_stop(1)
       t = timer_read(1)

#---------------------------------------------------------------------
#      perform verification and print results
#---------------------------------------------------------------------

       verify(dt, ss, sr, b_size, proc_num_zones, proc_zone_id, rho_i, us, vs, ws, speed, qs, square, rhs, forcing, u, nx, ny, nz)

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
