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

const tsum = Array{FloatType}(undef, t_last+2)
const tming = Array{FloatType}(undef, t_last+2)
const tmaxg = Array{FloatType}(undef, t_last+2)
const trecs = Array{FloatType}(undef, t_last+2)

const npbversion="3.4.2"

const t_names = ["total", "rhsx", "rhsy", "rhsz", "rhs", "xsolve", "ysolve", 
                 "zsolve", "txinvr", "pinvr", "ninvr", "tzetar", "add", 
                 "qbc_send", "qbc_recv", "bpack" , "exch" , "xcomm", "ycomm", 
                 "zcomm", "last"]
                 
# "go" is called "everywhere" in a distributed environment, where each worker is deployed in the master node of a cluster

function go(class::CLASS; itimer=0, npb_verbose=0, zone_mapping = nothing)
   
   niter = sp_class[class].niter
   dt    = sp_class[class].dt
   ratio = sp_class[class].ratio

   x_zones = sp_class[class].x_zones
   y_zones = sp_class[class].y_zones
   gx_size = sp_class[class].gx_size
   gy_size = sp_class[class].gy_size
   gz_size = sp_class[class].gz_size

   problem_size = sp_class[class].problem_size


   # =====================================================
   # the commands below are supposed to be prior executed  
   # @everywhere workers() using NPBApps              
   # @everywhere workers() using MPIClusterManagers   
   # =====================================================
   #for (w,np) in process_count
   #   fetch(@spawnat w addprocs(MPIWorkerManager(np)))
   #end

   clusters = Array{Tuple{Int64,Int64}}(undef, nworkers(role=:master)) 
   clusters .= map(w -> (w, @fetchfrom w nworkers(role=:master)), workers(role=:master))
  
   global num_clusters = length(clusters) 
   global max_zones = x_zones*y_zones
   num_zones = max_zones
   global max_threads = 128

   alloc_zone_vectors(x_zones, y_zones)

   zone_setup(x_zones, y_zones, gx_size, gy_size, gz_size, nx, nxmax, ny, nz, ratio, npb_verbose) 

   alloc_proc_space(num_clusters, max_zones) 

   num_processes .= map(c -> c[2], clusters)
   total_processes = reduce(+, num_processes)

   if isnothing(zone_mapping)
      proc_num_zones, proc_zone_id, zone_proc_id = map_zones(num_clusters, x_zones, y_zones, num_zones, nx, ny, nz, 1, total_processes, npb_verbose)    
   else
      proc_num_zones, proc_zone_id, zone_proc_id = map_zones(num_clusters, num_zones, zone_mapping) 
   end

 
   @printf(stdout, "\n\n NAS Parallel Benchmarks 3.4 -- SP Benchmark\n\n", )
   @printf(stdout, " Class %s\n", class)
   @printf(stdout, " Number of clusters: %4i\n", num_clusters)
   for clusterid = 1:num_clusters
      @printf(stdout, "  Cluster %4i with %4i processes\n", clusterid, num_processes[clusterid])
      for iz in 1:proc_num_zones[clusterid]
         zone_no = proc_zone_id[clusterid][iz]
         @printf(stdout, "  Zone %4i - Size: %4ix%4ix%4i \n", zone_no, nx[iz], ny[iz], nz[iz])
      end
   end
   @printf(stdout, " Iterations: %4i    dt: %11.7F\n", niter, dt)
   @printf(stdout, " Total number of processes: %6i\n", sum(num_processes))
   println(stdout) 

   global tmax_all = DataFlowVector(Array{FloatType}(undef, num_clusters))
   global timers_all = DataFlowVector(Array{Tuple{Vector{FloatType}, Vector{FloatType}, Vector{FloatType}}}(undef, num_clusters))

   global xce = zeros(5)
   global xcr = zeros(5)

   verified = false
   t = Threads.@spawn verified = verify(class, dt, niter)

   @everywhere workers() SP.go_cluster($clusters, $niter, $dt, $ratio, $x_zones, $y_zones, $gx_size, $gy_size, $gz_size, $problem_size, $nxmax, $nx, $ny, $nz, $proc_num_zones, $x_size, $y_size, $zone_proc_id, $proc_zone_id, $iz_west, $iz_east, $iz_south, $iz_north, $itimer, $npb_verbose)       

   wait(t)

   tmax = get(tmax_all, 1)
   tsum .= get(timers_all,1)[1] 
   tming .= get(timers_all,1)[2]
   tmaxg .= get(timers_all,1)[3]
   for i = 2:num_clusters
       tmax = max(get(tmax_all, i), tmax)
       tsum .= get(timers_all,1)[1] .+ tsum
       tmaxg .= max.(get(timers_all,1)[2], tmaxg)
       tming .= min.(get(timers_all,1)[3], tming)
   end

   mflops = 0.0e0
   if tmax != 0.
      for zone = 1:num_zones
        n3 = float(nx[zone])*ny[zone]*nz[zone]
        navg = (nx[zone] + ny[zone] + nz[zone])/3.0
        nsur = (nx[zone]*ny[zone] + nx[zone]*nz[zone] +
                ny[zone]*nz[zone])/3.0
        mflops = mflops + Float32( niter ) * 1.0e-6 *(881.174e0 * n3 - 4683.91e0 * nsur + 11484.5e0 * navg - 19272.4e0)/tmax
      end
   end

   print_results("SP-MZ", class, gx_size, gy_size, gz_size,
                     niter, num_clusters, num_processes, tmax, mflops, 
                     "          floating point",
                     verified, npbversion)

#---------------------------------------------------------------------
#      More timers
#---------------------------------------------------------------------

   @printf(stdout, " nprocs =%6i           minimum     maximum     average\n", num_clusters)
   for i = 1:t_last
      tsum[i] = tsum[i] / num_clusters
      @printf(stdout, " timer %2i(%8s) :  %10.4F  %10.4F  %10.4F\n", i, t_names[i], tming[i], tmaxg[i], tsum[i])
   end

   if itimer >= 2
      for ip = 0:num_clusters-1
         trecs = get(timers_all, ip+1)[1] ./ clusters[ip + 1][2]
         if (tmax == 0.0) tmax = 1.0 end
         @printf(stdout, "\n clusterid =%6i num_processes =%4i\n  SECTION   Time (secs)\n", ip, clusters[ip + 1][2])
         for i = 1:t_last
            @printf(stdout, "  %8s:%9.3F  (%6.2F)\n", t_names[i], trecs[i], trecs[i]*100.0/tmax)
            if i == t_rhs
               t = trecs[t_rhsx] + trecs[t_rhsy] + trecs[t_rhsz]
               @printf(stdout, "    --> total %8s:%9.3F  (%6.2F)\n", "sub-rhs", t, t*100.0/tmax)
               t = trecs[t_rhs] - t
               @printf(stdout, "    --> total %8s:%9.3F  (%6.2F)\n", "rest-rhs", t, t*100.0/tmax)
            elseif i == t_rdis2
               t = trecs[t_rdis1] + trecs[t_rdis2]
               @printf(stdout, "    --> total %8s:%9.3F  (%6.2F)\n", "exch_qbc", t, t*100.0/tmax)
            end
         end
      end
   end
end



function reportMaxTimeCluster(clusterid, tmax)
   set(tmax_all, clusterid+1, tmax)
end

function reportVerifiedCluster(clusterid, verified)
   set(verified_all, clusterid+1, verified)
end

function reportTimersCluster(clusterid, tsum, tming, tmaxg)
   set(timers_all, clusterid+1, (tsum, tming, tmaxg))
end

function send_face_through_driver(target_id, target_zone, face_data, _::Val{1})
   remotecall(send_east_face, target_id, target_zone, face_data; role = :master)
end

function send_face_through_driver(target_id, target_zone, face_data, _::Val{2})
   remotecall(send_west_face, target_id, target_zone, face_data; role = :master)
end

function send_face_through_driver(target_id, target_zone, face_data, _::Val{3})
   remotecall(send_north_face, target_id, target_zone, face_data; role = :master)
end

function send_face_through_driver(target_id, target_zone, face_data, _::Val{4})
   remotecall(send_south_face, target_id, target_zone, face_data; role = :master)
end