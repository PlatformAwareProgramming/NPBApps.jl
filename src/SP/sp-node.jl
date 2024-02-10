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



const t_names = ["total", "rhsx", "rhsy", "rhsz", "rhs", "xsolve", "ysolve", 
                 "zsolve", "txinvr", "pinvr", "ninvr", "tzetar", "add", 
                 "qbc_copy", "qbc_comm", "bpack" , "exch" , "xcomm", "ycomm", 
                 "zcomm", "last"]

function perform(clusterid_, clusters, niter, dt, ratio, x_zones, y_zones, gx_size, gy_size, gz_size, nxmax, nx, ny, nz, proc_num_zones, itimer_=false, npb_verbose_=0)

       global clusterid = clusterid_
       global num_clusters = length(clusters) 
       global itimer = itimer_
       global timeron = itimer > 0
       global npb_verbose = npb_verbose_
       global max_zones = x_zones * y_zones
       num_zones = max_zones
       
       setup_mpi(proc_num_zones) 

       if (!active) @goto L999 end

       class = set_class(niter, x_zones, y_zones, gx_size, gy_size, gz_size)

       alloc_field_space_zones(proc_num_zones)
       
       for zone = 1:proc_num_zones

          alloc_field_space(zone, [nx[zone], ny[zone], nz[zone]])
          make_set(zone, [nx[zone], ny[zone], nz[zone]])

          for c = 1:ncells
             if (cell_size[zone][1, c] > IMAX[zone]) ||(cell_size[zone][2, c] > JMAX[zone]) ||(cell_size[zone][3, c] > KMAX[zone])
                println(stdout, node, c, view(cell_size[zone], 1:3, c)...)
                println(stdout, " Problem size too big for compiled array sizes")
                @goto L999
             end
         end
      end

       for i = 1:t_last
          timer_clear(i)
       end

       for zone = 1:proc_num_zones         
         set_constants(zone, dt, gx_size, gy_size, gz_size, x_zones, y_zones) 
         initialize(zone) 
         lhsinit(zone) 
         exact_rhs(zone) 
         compute_buffer_size(zone, 5)
      end


      if (no_nodes > 1)
         ss = SA[start_send_east::Int start_send_west::Int start_send_north::Int start_send_south::Int start_send_top::Int start_send_bottom::Int]
         sr = SA[start_recv_east::Int start_recv_west::Int start_recv_north::Int start_recv_south::Int start_recv_top::Int start_recv_bottom::Int]
         b_size = SA[east_size::Int west_size::Int north_size::Int south_size::Int top_size::Int bottom_size::Int]
      else
         ss = nothing
         sr = nothing
         b_size = nothing
      end

       requests = Array{Array{MPI.Request}}(undef,proc_num_zones)
       s = Array{Array{Float64}}(undef,proc_num_zones)
       for iz = 1:proc_num_zones
         requests[iz] = Array{MPI.Request}(undef,12)
         s[iz] = Array{Float64}(undef,5)
       end

#---------------------------------------------------------------------
#      do one time step to touch all code, and reinitialize
#---------------------------------------------------------------------
       for zone = 1:proc_num_zones
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
               ss,
               sr,
               b_size,
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

       for STEP = 1:niter

          if node == root
             if mod(STEP, 20) == 0 || STEP == 1
                @printf(stdout, "%2i: Time step %4i\n", clusterid,   STEP)
              end
          end

          for zone = 1:proc_num_zones
                          
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
                  ss,
                  sr,
                  b_size,
                  s[zone],
                  requests[zone]
            )
         end
       end

       timer_stop(1)
       t = timer_read(1)

#---------------------------------------------------------------------
#      perform verification and print results
#---------------------------------------------------------------------

       verify(dt, ss, sr, b_size, proc_num_zones, rho_i, us, vs, ws, speed, qs, square, rhs, forcing, u, nx, ny, nz)

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
