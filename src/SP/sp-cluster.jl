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


function update_and_test_face_out_count(f, z, l1, h1, l2, h2, n1, n2)
   t = false
   lock(lkout)
   try   
         face_out_count[z][f] += (h1 - l1 + 1)*(h2 - l2 + 1)
         t = face_out_count[z][f] == (n1[z]-2)*(n2[z]-2)*5 
         if t
            face_out_count[z][f] = 0
         end
   finally
         unlock(lkout)
   end
   return t
end

function deposit_face(z, m, l1, h1, l2, h2, buffer, _::Val{1})
   view(faces_out[z][1], l1:h1, l2:h2, m) .= buffer 
   lock(lk1)
   try
      if update_and_test_face_out_count(1, z, l1, h1, l2, h2, ny, nz) 
         westiz = iz_west[proc_zone_id[clusterid + 1][z]]
         westid = zone_proc_id[westiz] 
         if westid != clusterid  
            remotecall(send_east_face, westid + 2, z, faces_out[z][1]; role=:worker)
         else
            send_east_face(z, faces_out[z][1])
         end
      end
   finally
      unlock(lk1)
   end
end

function deposit_face(z, m, l1, h1, l2, h2, buffer, _::Val{2})
   view(faces_out[z][2], l1:h1, l2:h2, m) .= buffer 
   lock(lk1)
   try
      if update_and_test_face_out_count(2, z, l1, h1, l2, h2, ny, nz)
         eastiz = iz_east[proc_zone_id[clusterid + 1][z]]
         eastid = zone_proc_id[eastiz] 
         if eastid != clusterid
            remotecall(send_west_face, eastid + 2, z, faces_out[z][2]; role=:worker)
         else
            send_west_face(z, faces_out[z][2])
         end
      end
   finally
      unlock(lk1)
   end
end

function deposit_face(z, m, l1, h1, l2, h2, buffer, _::Val{3})
   view(faces_out[z][3], l1:h1, l2:h2, m) .= buffer 
   lock(lk1)
   try
      if update_and_test_face_out_count(3, z, l1, h1, l2, h2, nx, nz)
         southiz = iz_south[proc_zone_id[clusterid + 1][z]]
         southid = zone_proc_id[southiz]
         if southid != clusterid 
            remotecall(send_north_face, southid + 2, z, faces_out[z][3]; role=:worker)
         else 
            send_north_face(z, faces_out[z][3])
         end
      end
   finally
      unlock(lk1)
   end
end

function deposit_face(z, m, l1, h1, l2, h2, buffer, _::Val{4})
   view(faces_out[z][4], l1:h1, l2:h2, m) .= buffer 
   lock(lk1)
   try
      if update_and_test_face_out_count(4, z, l1, h1, l2, h2, nx, nz)
         northiz = iz_north[proc_zone_id[clusterid + 1][z]]
         northid = zone_proc_id[northiz]
         if northid != clusterid  
            remotecall(send_south_face, northid + 2, z, faces_out[z][4]; role=:worker)
         else
            send_south_face(z, faces_out[z][4])
         end
      end
   finally
      unlock(lk1)
   end
end


function send_west_face(z, buffer)
   faces_in[z][1] .= buffer 
   lock(faces_in_cond[z][1])
   try
      face_in_count[z][1] = 0
      notify(faces_in_cond[z][1])
   finally
      unlock(faces_in_cond[z][1])
   end
end

function send_east_face(z, buffer)
   faces_in[z][2] .= buffer 
   lock(faces_in_cond[z][2])
   try
      face_in_count[z][2] = 0
      notify(faces_in_cond[z][2])
   finally
      unlock(faces_in_cond[z][2])
   end
end

function send_north_face(z, buffer)
   faces_in[z][4] .= buffer 
   lock(faces_in_cond[z][4])
   try
      face_in_count[z][4] = 0
      notify(faces_in_cond[z][4])
   finally
      unlock(faces_in_cond[z][4])
   end
end

function send_south_face(z, buffer)
   faces_in[z][3] .= buffer 
   lock(faces_in_cond[z][3])
   try
      face_in_count[z][3] = 0
      notify(faces_in_cond[z][3])
   finally
      unlock(faces_in_cond[z][3])
   end
end

function update_face_in_count(z, l1, h1, l2, h2, f)
   lock(lkin)
   try
      face_in_count[z][f] += (h1-l1+1)*(h2-l2+1)
   finally
      unlock(lkin)
   end
end

function collect_face(z, m, l1, h1, l2, h2, _::Val{1})
   lock(lk2)
   try
      lock(faces_in_cond[z][1])
      try
         while face_in_count[z][1] == (ny[z]-2)*(nz[z]-2)*5
            wait(faces_in_cond[z][1])
         end
      finally
         unlock(faces_in_cond[z][1])
      end
      r = view(faces_in[z][1], l1:h1, l2:h2, m)    
      update_face_in_count(z, l1, h1, l2, h2, 1)
      r
   finally
      unlock(lk2)
   end
end

function collect_face(z, m, l1, h1, l2, h2, _::Val{2})
   lock(lk2)
   try
      lock(faces_in_cond[z][2])
      try
         while face_in_count[z][2] == (ny[z]-2)*(nz[z]-2)*5
            wait(faces_in_cond[z][2])
         end
      finally
         unlock(faces_in_cond[z][2])
      end
      r = view(faces_in[z][2], l1:h1, l2:h2, m) 
      update_face_in_count(z, l1, h1, l2, h2, 2)
      r
   finally
      unlock(lk2)
   end
end

function collect_face(z, m, l1, h1, l2, h2, _::Val{3})
   lock(lk2)
   try
      lock(faces_in_cond[z][3])
      try
         while face_in_count[z][3] == (nx[z]-2)*(nz[z]-2)*5 
            wait(faces_in_cond[z][3])
         end
      finally
         unlock(faces_in_cond[z][3])
      end
      r = view(faces_in[z][3], l1:h1, l2:h2, m) 
      update_face_in_count(z, l1, h1, l2, h2, 3)
      r
   finally
      unlock(lk2)
   end
end

function collect_face(z, m, l1, h1, l2, h2, _::Val{4})
   lock(lk2)
   try
      lock(faces_in_cond[z][4])
      try
         while face_in_count[z][4] == (nx[z]-2)*(nz[z]-2)*5 
            wait(faces_in_cond[z][4])
         end
      finally
         unlock(faces_in_cond[z][4])
      end
      r = view(faces_in[z][4], l1:h1, l2:h2, m) 
      update_face_in_count(z, l1, h1, l2, h2, 4)
      r
   finally
      unlock(lk2)
   end
end

function reportMaxTimeNode(tmax)
   remotecall(reportMaxTimeCluster, 1, clusterid, tmax; role = :worker)
end

function reportVerifiedNode(verified)
   remotecall(reportVerifiedCluster, 1, clusterid, verified; role = :worker)
end

function reportTimersNode(tsum, tming, tmaxg)
   remotecall(reportTimersCluster, 1, clusterid, tsum, tming, tmaxg; role = :worker)
end

function go_cluster(clusters, niter, dt, ratio, x_zones, y_zones, gx_size, gy_size, gz_size, nxmax, nx_, ny_, nz_, proc_num_zones_all, x_size, y_size, zone_proc_id_, proc_zone_id_, iz_west_, iz_east_, iz_south_, iz_north_, itimer, npb_verbose)
   global clusterid = Distributed.myid(role=:worker) - 2 

   global nx = nx_
   global ny = ny_
   global nz = nz_

   global zone_proc_id = zone_proc_id_
   global proc_zone_id = proc_zone_id_
   global iz_west = iz_west_ 
   global iz_east = iz_east_ 
   global iz_south = iz_south_
   global iz_north = iz_north_

   proc_num_zones = proc_num_zones_all[clusterid+1]

   global faces_out = Array{Array{Array{Float64}}}(undef, proc_num_zones)
   for iz = 1:proc_num_zones
      faces_out[iz] = Array{Array{Float64}}(undef, 4)
      faces_out[iz][1] = Array{Float64}(undef, ny[iz], gz_size, 5) 
      faces_out[iz][2] = Array{Float64}(undef, ny[iz], gz_size, 5) 
      faces_out[iz][3] = Array{Float64}(undef, nx[iz], gz_size, 5) 
      faces_out[iz][4] = Array{Float64}(undef, nx[iz], gz_size, 5) 
   end

   global faces_in = Array{Array{Array{Float64}}}(undef, proc_num_zones)
   for iz = 1:proc_num_zones
      faces_in[iz] = Array{Array{Float64}}(undef, 4)
      faces_in[iz][1] = Array{Float64}(undef, ny[iz], gz_size, 5)
      faces_in[iz][2] = Array{Float64}(undef, ny[iz], gz_size, 5)
      faces_in[iz][3] = Array{Float64}(undef, nx[iz], gz_size, 5)
      faces_in[iz][4] = Array{Float64}(undef, nx[iz], gz_size, 5)
   end

   global face_in_count = Array{Array{Int64}}(undef, proc_num_zones)
   global face_out_count = Array{Array{Int64}}(undef, proc_num_zones)
   global faces_in_cond = Array{Array{Threads.Condition}}(undef, proc_num_zones)
   for iz = 1:proc_num_zones
      face_in_count[iz] = Array{Int64}(undef, 6)
      face_out_count[iz] = zeros(Int64, 6)
      faces_in_cond[iz] = Array{Threads.Condition}(undef, 6)
      for f = 1:4
         faces_in_cond[iz][f] = Threads.Condition()
      end
      face_in_count[iz][1] = (ny[iz]-2)*(nz[iz]-2)*5
      face_in_count[iz][2] = (ny[iz]-2)*(nz[iz]-2)*5
      face_in_count[iz][3] = (nx[iz]-2)*(nz[iz]-2)*5 
      face_in_count[iz][4] = (nx[iz]-2)*(nz[iz]-2)*5
   end

   global lkout = ReentrantLock()
   global lkin = ReentrantLock()

   global lk1 = ReentrantLock()
   global lk2 = ReentrantLock()

   SP.go_node(clusterid, clusters, niter, dt, ratio, x_zones, y_zones, gx_size, gy_size, gz_size, nxmax, nx, ny, nz, proc_num_zones, itimer, npb_verbose) 
end

function go_node(cluster_id, no_nodes, niter, dt, ratio, x_zones, y_zones, gx_size, gy_size, gz_size, nxmax, nx, ny, nz, proc_num_zones, itimer, npb_verbose)
   @everywhere workers() SP.perform($cluster_id, $no_nodes, $niter, $dt, $ratio, $x_zones, $y_zones, $gx_size, $gy_size, $gz_size, $nxmax, $nx, $ny, $nz, $proc_num_zones, $itimer, $npb_verbose) 
end


