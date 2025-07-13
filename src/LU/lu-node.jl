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



function perform(clusterid_, clusters, itmax, inorm, dt, ratio, x_zones, y_zones, gx_size, gy_size, gz_size, problem_size, nxmax, nx0, ny0, nz0, proc_num_zones, proc_zone_id_, zone_proc_id_, iz_west, iz_east, iz_south, iz_north, itimer_=false, npb_verbose_=0)

      global clusterid = clusterid_
      global num_clusters = length(clusters) 
      global itimer = itimer_
      global timeron = itimer > 0
      global npb_verbose = npb_verbose_
      global max_zones = x_zones * y_zones
      global proc_zone_id = proc_zone_id_
      global zone_proc_id = zone_proc_id_
      num_zones = max_zones
      niter = itmax
      
      omega = omega_default
      tolrsd = tolrsddef

#---------------------------------------------------------------------
#   initialize communications
#---------------------------------------------------------------------
      init_comm(proc_num_zones) 

      if (!active) @goto L999 end

      class = set_class(itmax, gx_size, gy_size, gz_size)

      for i = 1:t_last
         timer_clear(i)
      end

#---------------------------------------------------------------------
#   allocate space
#---------------------------------------------------------------------
      alloc_field_space_zones(proc_num_zones)

      nx = Array{Int64}(undef, proc_num_zones)
      ny = Array{Int64}(undef, proc_num_zones)
      nz = Array{Int64}(undef, proc_num_zones)
      ipt = Array{Int64}(undef, proc_num_zones)
      jpt = Array{Int64}(undef, proc_num_zones)
      ist = Array{Int64}(undef, proc_num_zones)
      jst = Array{Int64}(undef, proc_num_zones)
      iend = Array{Int64}(undef, proc_num_zones)
      jend = Array{Int64}(undef, proc_num_zones)

      row = Array{Int64}(undef, proc_num_zones)
      col = Array{Int64}(undef, proc_num_zones)

      south = Array{Int64}(undef, proc_num_zones)
      east = Array{Int64}(undef, proc_num_zones)
      north = Array{Int64}(undef, proc_num_zones)
      west = Array{Int64}(undef, proc_num_zones)

      south2 = Array{Int64}(undef, proc_num_zones)
      east2 = Array{Int64}(undef, proc_num_zones)
      north2 = Array{Int64}(undef, proc_num_zones)
      west2 = Array{Int64}(undef, proc_num_zones)

      for iz = 1:proc_num_zones
           
            zone = proc_zone_id[iz]

#---------------------------------------------------------------------
#   set up processor grid
#---------------------------------------------------------------------
            row[iz], col[iz] = proc_grid(nx0[zone], ny0[zone], nz0[zone])

#---------------------------------------------------------------------
#   determine the neighbors
#---------------------------------------------------------------------
            south[iz], east[iz], north[iz], west[iz], south2[iz], east2[iz], north2[iz], west2[iz] = neighbors(row[iz], col[iz]) 

            @info "$clusterid/$node: zone=$zone --- westzone=$(iz_west[zone]) eastzone=$(iz_east[zone]) northzone=$(iz_north[zone]) southzone=$(iz_south[zone])"
            @info "$clusterid/$node: zone=$zone --- west=$(west[iz]) east=$(east[iz]) north=$(north[iz]) south=$(south[iz])"
            @info "$clusterid/$node: zone=$zone --- west2=$(west2[iz]) east2=$(east2[iz]) north2=$(north2[iz]) south2=$(south2[iz])"
#---------------------------------------------------------------------
#   set up sub-domain sizes (calculate nx, ny, nz, ipt, jpt, ist, jst, iend, jend for zone iz)
#---------------------------------------------------------------------
            subdomain(iz, row[iz], col[iz], west[iz], east[iz], north[iz], south[iz], nx0[zone], ny0[zone], nz0[zone], nx, ny, nz, ipt, jpt, ist, jst, iend, jend) 

            @info "$clusterid/$node: zone=$zone --- nx=$(nx[iz]), ny=$(ny[iz]), nz=$(nz[iz]), ipt=$(ipt[iz]), jpt=$(jpt[iz]), ist=$(ist[iz]), jst=$(jst[iz]), iend=$(iend[iz]), jend=$(jend[iz])"

            alloc_field_space(iz, #=nx0[zone], ny0[zone], nz0[zone],=# nx[iz], ny[iz], nz[iz], problem_size)

#---------------------------------------------------------------------
#   set up coefficients
#---------------------------------------------------------------------
            setcoeff(nx0[zone], ny0[zone], nz0[zone])

#---------------------------------------------------------------------
#   set the boundary values for dependent variables
#---------------------------------------------------------------------
            setbv(u[iz], nx0[zone], ny0[zone], nz0[zone], nx[iz], ny[iz], nz[iz], west[iz], east[iz], north[iz], south[iz], ipt[iz], jpt[iz])  

#---------------------------------------------------------------------
#   set the initial values for dependent variables
#---------------------------------------------------------a------------
            setiv(u[iz], nx0[zone], ny0[zone], nz0[zone], nx[iz], ny[iz], nz[iz], ipt[iz], jpt[iz])  

#---------------------------------------------------------------------
#   compute the forcing term based on prescribed exact solution
#---------------------------------------------------------------------
            erhs(frct[iz], rsd[iz], flux[iz], buf[iz], buf1[iz], nx0[zone], ny0[zone], nz0[zone], nx[iz], ny[iz], nz[iz], west[iz], east[iz], north[iz], south[iz], ipt[iz], jpt[iz], ist[iz], jst[iz], iend[iz], jend[iz], comm_solve[iz])

#---------------------------------------------------------------------
#   initialize a,b,c,d to zero (guarantees that page tables have been
#   formed, if applicable on given architecture, before timestepping).
#---------------------------------------------------------------------
            init_workarray(nx[iz], a[iz], b[iz], c[iz], d[iz])

#---------------------------------------------------------------------
#   compute the steady-state residuals
#---------------------------------------------------------------------
            rhs(               
                  comm_solve[iz], 
                  u[iz],
                  rsd[iz],
                  frct[iz],
                  flux[iz],
                  buf[iz],
                  buf1[iz],
                  south[iz],
                  east[iz],
                  north[iz],
                  west[iz],
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
                  nx[iz],
                  ny[iz],
                  nz[iz],
                  ist[iz],
                  iend[iz],
                  jst[iz],
                  jend[iz],
               )

      end


     requests = Array{MPI.Request}(undef,8*proc_num_zones)

     # if no_nodes > 1 && (num_clusters > 1 || proc_num_zones > 1)
         exch_qbc(proc_num_zones, zone_proc_id, proc_zone_id, u, row, col, west, east, north, south, west2, east2, north2, south2, nx, ny, nz, timeron, ist, iend, jst, jend, iz_west, iz_east, iz_south, iz_north, comm_exch, buf_exch_w_out, buf_exch_e_out, buf_exch_n_out, buf_exch_s_out, buf_exch_w_in, buf_exch_e_in, buf_exch_n_in, buf_exch_s_in, requests)
     # end


     #@info "- - - - - before"
     #for iz = 1:proc_num_zones
     #    zone = proc_zone_id[iz]
     #   write_u(iz, u[iz], nx0[zone], ny0[zone], nz0[zone], nx[iz], ny[iz], nz[iz], ipt[iz], jpt[iz]) 
     #end

#---------------------------------------------------------------------
#   perform one SSOR iteration to touch all data and program pages 
#---------------------------------------------------------------------
      Threads.@threads for iz = 1:proc_num_zones
            zone = proc_zone_id[iz]
            
            ssor( iz, ipt[iz], jpt[iz], proc_num_zones,
                  comm_solve[iz],
                  u[iz],
                  rsd[iz],
                  frct[iz],
                  flux[iz],
                  a[iz],
                  b[iz],
                  c[iz],
                  d[iz],
                  buf[iz], 
                  buf1[iz],
                  south[iz],
                  east[iz],
                  north[iz],
                  west[iz],
                  dt,
                  omega,
                  nx0[zone],
                  ny0[zone],
                  nz0[zone],
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
                  nx[iz],
                  ny[iz],
                  nz[iz],
                  ist[iz],
                  iend[iz],
                  jst[iz],
                  jend[iz],
            )
      end

      #@info "- - - - - after"
      #for iz = 1:proc_num_zones
      #    zone = proc_zone_id[iz]
      #    write_u(iz, u[iz], nx0[zone], ny0[zone], nz0[zone], nx[iz], ny[iz], nz[iz], ipt[iz], jpt[iz]) 
      #end


      #---------------------------------------------------------------------
      #   reset the boundary and initial values
      #---------------------------------------------------------------------
      Threads.@threads for iz = 1:proc_num_zones
            zone = proc_zone_id[iz]

            setbv(u[iz], nx0[zone], ny0[zone], nz0[zone], nx[iz], ny[iz], nz[iz], west[iz], east[iz], north[iz], south[iz], ipt[iz], jpt[iz])
            setiv(u[iz], nx0[zone], ny0[zone], nz0[zone], nx[iz], ny[iz], nz[iz], ipt[iz], jpt[iz])

            #---------------------------------------------------------------------
            #   initialize a,b,c,d to zero (guarantees that page tables have been
            #   formed, if applicable on given architecture, before timestepping).
            #---------------------------------------------------------------------
            init_workarray(nx[iz], a[iz], b[iz], c[iz], d[iz])

            #---------------------------------------------------------------------
            #   compute the steady-state residuals
            #---------------------------------------------------------------------
            rhs(               
                  comm_solve[iz], 
                  u[iz],
                  rsd[iz],
                  frct[iz],
                  flux[iz],
                  buf[iz],
                  buf1[iz],
                  south[iz],
                  east[iz],
                  north[iz],
                  west[iz],
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
                  nx[iz],
                  ny[iz],
                  nz[iz],
                  ist[iz],
                  iend[iz],
                  jst[iz],
                  jend[iz],
               )

       end

       for i = 1:t_last
           timer_clear(i)
       end
   
       MPI.Barrier(comm_setup)
   
       timer_clear(1)
       timer_start(1)
   
       Q = 1

       timer_clear(64); t_64 = 0.0; t_64s = 0.0
       timer_clear(63); t_63 = 0.0; t_63s = 0.0

      # @info "$clusterid/$node: GO! NUM_THREADS = $(Threads.nthreads())"

      #---------------------------------------------------------------------
      #   perform the SSOR iterations
      #---------------------------------------------------------------------
      for istep = 1:niter

            if node == root
               Q = istep > 1 && t_64 < 5.0 ? ceil(5.0 / t_64) : Q

               if mod(istep, Q) == 0 || istep == itmax || istep == 1
                  if (niter > 1) 
                     @printf(stdout, "%2i: Time step %4i  -- %12.5F  -- %12.5F --- %4i \n", clusterid, istep, t_63s/(istep-1), t_64s/(istep-1), Q)
                  end
               end
            end

            if istep == 30
                  # RESTART THE TIMERS TO MINIMIZE TRANSIENT EFFECTS
                  @info "restarting timers"
                  timer_clear(64); t_64 = 0.0; t_64s = 0.0
                  timer_clear(63); t_63 = 0.0; t_63s = 0.0
            end

           timer_start(63)
         
           # if no_nodes > 1 && (num_clusters > 1 || proc_num_zones > 1)
               exch_qbc(proc_num_zones, zone_proc_id, proc_zone_id, u, row, col, west, east, north, south, west2, east2, north2, south2, nx, ny, nz, timeron, ist, iend, jst, jend, iz_west, iz_east, iz_south, iz_north, comm_exch, buf_exch_w_out, buf_exch_e_out, buf_exch_n_out, buf_exch_s_out, buf_exch_w_in, buf_exch_e_in, buf_exch_n_in, buf_exch_s_in, requests) 
           # end

            t_63 = timer_stop(63); t_63s += t_63

            timer_start(64)

         #   if istep == 2
         #      @info "- - - - - aaa $istep"
         #      for iz = 1:proc_num_zones
         #         zone = proc_zone_id[iz]
         #         write_u(iz, u[iz], nx0[zone], ny0[zone], nz0[zone], nx[iz], ny[iz], nz[iz], ipt[iz], jpt[iz]) 
         #      end
         #   end

            Threads.@threads for iz = 1:proc_num_zones
                  zone = proc_zone_id[iz]
                  ssor( iz, ipt[iz], jpt[iz], proc_num_zones,
                        comm_solve[iz],
                        u[iz],
                        rsd[iz],
                        frct[iz], 
                        flux[iz],
                        a[iz],
                        b[iz], 
                        c[iz],
                        d[iz],
                        buf[iz],
                        buf1[iz],
                        south[iz],
                        east[iz],
                        north[iz],
                        west[iz],
                        dt,
                        omega,
                        nx0[zone],
                        ny0[zone],
                        nz0[zone],
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
                        nx[iz],
                        ny[iz],
                        nz[iz],
                        ist[iz],
                        iend[iz],
                        jst[iz],
                        jend[iz],
                        )
                  end

         
            #if istep == 2
            #   @goto L999
            #end

            t_64 = timer_stop(64); t_64s += t_64
     end

      timer_stop(1)
      wtime = timer_read(1)      

      errnm_all = zeros(5)
      rsdnm_all = zeros(5)

      frc = 0.0
      Threads.@threads for iz = 1:proc_num_zones
          zone = proc_zone_id[iz]

          #---------------------------------------------------------------------
          #   compute the solution error
          #---------------------------------------------------------------------
          ERROR(u[iz], nx0[zone], ny0[zone], nz0[zone], ipt[iz], jpt[iz], ist[iz], jst[iz], iend[iz], jend[iz], errnm[iz], comm_solve[iz])

          #---------------------------------------------------------------------
          #   compute the max-norms of newton iteration residuals
          #---------------------------------------------------------------------
          l2norm( nx0[zone], ny0[zone], nz0[zone], ist[iz], iend[iz], jst[iz], jend[iz], rsd[iz], rsdnm[iz], comm_solve[iz], timeron,)

          errnm_all .+= errnm[iz] # xce
          rsdnm_all .+= rsdnm[iz] # xcr

          #---------------------------------------------------------------------
          #   compute the surface integral
          #---------------------------------------------------------------------
          frc0 =  pintgr(u[iz], phi1[iz], phi2[iz], ny0[zone], nz0[zone], nx[iz], ny[iz], nz[iz], ipt[iz], jpt[iz], west[iz], east[iz], south[iz], north[iz], comm_solve[iz])       
          frc += frc0

      end

#---------------------------------------------------------------------
#   verification test
#---------------------------------------------------------------------
      verify(rsdnm_all, errnm_all, frc)

      tmax = MPI.Reduce(wtime, MPI.MAX, root, comm_setup)
      if node == root
         remotecall(reportMaxTimeNode, 1, tmax; role=:worker)
      end

#---------------------------------------------------------------------
#      More timers
#---------------------------------------------------------------------

      if (!timeron) @goto L999 end

      t1 = Array{FloatType}(undef, t_last+2)

      for i = 1:t_last
         t1[i] = timer_read(i)
      end
      t1[t_rhs] = t1[t_rhs] - t1[t_exch]
      t1[t_last+2] = t1[t_lcomm]+t1[t_ucomm]+t1[t_rcomm]+t1[t_exch]
      t1[t_last+1] = t1[t_total] - t1[t_last+2]

      tsum = MPI.Reduce(t1, MPI.SUM, 0, comm_setup)
      tming = MPI.Reduce(t1, MPI.MIN, 0, comm_setup)
      tmaxg = MPI.Reduce(t1, MPI.MAX, 0, comm_setup)

      if node == root
         tavg = tsum ./ no_nodes
         remotecall(reportTimersNode, 1, tavg, tming, tmaxg; role=:worker)
      end
   
      @label L999
      #MPI.Finalize()

      return nothing
end


function init_workarray(nx, a, b, c, d)
#---------------------------------------------------------------------
#   initialize a,b,c,d to zero (guarantees that page tables have been
#   formed, if applicable on given architecture, before timestepping).
#---------------------------------------------------------------------
for i = 1:nx
      for m = 1:5
         for k = 1:5
            a[k, m, i] = 0.0e0
            b[k, m, i] = 0.0e0
            c[k, m, i] = 0.0e0
            d[k, m, i] = 0.0e0
         end
      end
   end

end