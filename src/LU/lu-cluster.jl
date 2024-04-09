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


function update_and_test_face_out_count(f, z, l1, h1, l2, h2, n1, n2)
   t = false
   lock(lkout)
   try   
         face_out_count[z][f] += (h1 - l1 + 1)*(h2 - l2 + 1)*5
         t = face_out_count[z][f] == n1*n2*5 
         #@info "$clusterid: update_and_test_face_out_count f=$f t=$t face_out_count[$z][$f]=$(face_out_count[z][f]) --- $(n1*n2*5)"
         if t
            face_out_count[z][f] = 0
         end
   finally
         unlock(lkout)
   end
   return t
end

proc = ["send_east_face", "send_west_face", "send_north_face", "send_south_face"]

function perform_deposit_face(z, target_zone, l1, h1, l2, h2, f, n1, n2, buffer, send_proc)

   #@info "$clusterid: perform_deposit_face($f) BEGIN --- z=$z target_zone=$target_zone $l1/$h1/$l2/$h2"
   view(face_out[z][f], 1:5, l1:h1, l2:h2) .= buffer 
   if update_and_test_face_out_count(f, z, l1, h1, l2, h2, n1, n2)
      clusterid_ = zone_proc_id[target_zone] 
      if clusterid_ != clusterid  
         #@info "$clusterid: remotecall($(proc[f]), $(clusterid_ + 2), $target_zone, face_out[$z][$f]; role=:worker) - BEGIN"
         lock(lk1)
         try
            remotecall(send_proc, clusterid_ + 2, target_zone, face_out[z][f]; role=:worker)
         finally
            unlock(lk1)
         end
         #@info "$clusterid: remotecall($(proc[f]), $(clusterid_ + 2), $target_zone, face_out[$z][$f]; role=:worker) - END"
      else
         #@info "$clusterid: (local) $(proc[f])($target_zone, face_out[$z][$f]) - BEGIN"
         send_proc(target_zone, face_out[z][f])
         #@info "$clusterid: (local) $(proc[f])($target_zone, face_out[$z][$f]) - END"
      end
   end
   #@info "$clusterid: perform_deposit_face($f) END --- z=$z target_zone=$target_zone $l1/$h1/$l2/$h2"
end

# WEST
function deposit_face(z, l1, h1, l2, h2, buffer, _::Val{1})
   zone = proc_zone_id[z]
   perform_deposit_face(z, iz_west[proc_zone_id[z]], l1, h1, l2, h2, 1, nx[zone]#=-2=#, nz[zone]#=-2=#, buffer, send_east_face)
end

# EAST
function deposit_face(z, l1, h1, l2, h2, buffer, _::Val{2})
   zone = proc_zone_id[z]
   perform_deposit_face(z, iz_east[proc_zone_id[z]], l1, h1, l2, h2, 2, nx[zone]#=-2=#, nz[zone]#=-2=#, buffer, send_west_face)
end

# SOUTH
function deposit_face(z, l1, h1, l2, h2, buffer, _::Val{3})
   zone = proc_zone_id[z]
   perform_deposit_face(z, iz_south[proc_zone_id[z]], l1, h1, l2, h2, 3, ny[zone]#=-2=#, nz[zone]#=-2=#, buffer, send_north_face)
end

# NORTH
function deposit_face(z, l1, h1, l2, h2, buffer, _::Val{4})
   zone = proc_zone_id[z]
   perform_deposit_face(z, iz_north[proc_zone_id[z]], l1, h1, l2, h2, 4, ny[zone]#=-2=#, nz[zone]#=-2=#, buffer, send_south_face)
end

function proc_zone_id_inv(zone)
   for iz in 1:proc_num_zones
     proc_zone_id[iz] == zone && return iz
   end
   #@info "$clusterid: zone $zone not found"
end

function send_west_face(zone, buffer)
   z = proc_zone_id_inv(zone)
   #@info "$clusterid: send_west_face - BEGIN 1 --- zone=$zone z=$z"
   face_receive_wait(z, 1, nx[zone], nz[zone])
   #@info "$clusterid: send_west_face - BEGIN 2 --- zone=$zone z=$z"
   face_receive_and_notify(buffer, z, 1)
   #@info "$clusterid: send_west_face - END --- zone=$zone z=$z"
end

function send_east_face(zone, buffer)
   z = proc_zone_id_inv(zone)
   #@info "$clusterid: send_east_face - BEGIN 1 --- zone=$zone z=$z"
   face_receive_wait(z, 2, nx[zone], nz[zone])
   #@info "$clusterid: send_east_face - BEGIN 2 --- zone=$zone z=$z"
   face_receive_and_notify(buffer, z, 2)
   #@info "$clusterid: send_east_face - END --- zone=$zone z=$z "
end

function send_south_face(zone, buffer)
   z = proc_zone_id_inv(zone)
   #@info "$clusterid: send_south_face - BEGIN 1 --- zone=$zone z=$z"
   face_receive_wait(z, 3, ny[zone], nz[zone])
   #@info "$clusterid: send_south_face - BEGIN 2 --- zone=$zone z=$z"
   face_receive_and_notify(buffer, z, 3)
   #@info "$clusterid: send_south_face - END --- zone=$zone z=$z"
end

function send_north_face(zone, buffer)
   z = proc_zone_id_inv(zone)
   #@info "$clusterid: send_north_face - BEGIN 1 --- zone=$zone z=$z"
   face_receive_wait(z, 4, ny[zone], nz[zone])
   #@info "$clusterid: send_north_face - BEGIN 2 --- zone=$zone z=$z"
   face_receive_and_notify(buffer, z, 4)
   #@info "$clusterid: send_north_face - END --- zone=$zone z=$z"
end

function face_receive_wait(z, f, n1, n2)
   try
      lock(face_in_receive[z][f])
      while face_in_count[z][f] < (n1#=-2=#)*(n2#=-2=#)*5    
         #@info "$clusterid: receive wait BEGIN f=$f z=$z --- $(face_in_count[z][f])<$((n1-2)*(n2-2)*5)"   
         wait(face_in_receive[z][f])
         #@info "$clusterid: receive wait END f=$f z=$z --- $(face_in_count[z][f])<$((n1-2)*(n2-2)*5)"   
      end
   finally
      unlock(face_in_receive[z][f])
   end
end

function face_receive_and_notify(buffer, z, f)
   #@info "$clusterid: face_receive_and_notify - z=$z f=$f --- $buffer"
   face_in[z][f] .= buffer 
   lock(face_in_collect[z][f])
   try
      face_in_count[z][f] = 0
      notify(face_in_collect[z][f]#=; all=true=#)
   finally
      unlock(face_in_collect[z][f])
   end
end


function update_face_in_count(z, l1, h1, l2, h2, f, n1, n2)
   lock(lkin)
   try
      face_in_count[z][f] += (h1-l1+1)*(h2-l2+1)*5
   finally
      unlock(lkin)
   end

   lock(face_in_receive[z][f])
   try
      if face_in_count[z][f] == (n1#=-2=#)*(n2#=-2=#)*5
         #@info "$clusterid: face collect finished - notify ! z=$z f=$f count=$(face_in_count[z][f])"
         notify(face_in_receive[z][f]#=; all=true=#)
      end
   finally
      unlock(face_in_receive[z][f])
   end
end

function perform_collect_face(z, l1, h1, l2, h2, f, n1, n2)
   #@info "$clusterid: collect_face($f) BEGIN 1 --- z=$z $l1/$h1/$l2/$h2"
   #lock(lk2)
   #@info "$clusterid: collect_face($f) BEGIN 2 --- z=$z $l1/$h1/$l2/$h2"
   #try
      lock(face_in_collect[z][f])
      try
         while face_in_count[z][f] == n1*n2*5
            wait(face_in_collect[z][f])
         end
         #@info "$clusterid: collect_face($f) PASSED --- z=$z"
      finally
         unlock(face_in_collect[z][f])
      end
      
      r = face_in[z][f][1:5, l1:h1, l2:h2]
      @async update_face_in_count(z, l1, h1, l2, h2, f, n1, n2)
      #@info "$clusterid: collect_face($f) END --- z=$z -- $(face_in_count[z][f]) == $(n1*n2*5)"
      r
      #view(face_in[z][f], 1:5, l1:h1, l2:h2)   
   #finally
   #   unlock(lk2)
   #end
end

function collect_face(z, l1, h1, l2, h2, _::Val{1})
   zone = proc_zone_id[z]
   perform_collect_face(z, l1, h1, l2, h2, 1, nx[zone]#=-2=#, nz[zone]#=-2=#)
end

function collect_face(z, l1, h1, l2, h2, _::Val{2})
   zone = proc_zone_id[z]
   perform_collect_face(z, l1, h1, l2, h2, 2, nx[zone]#=-2=#, nz[zone]#=-2=#)
end

function collect_face(z, l1, h1, l2, h2, _::Val{3})
   zone = proc_zone_id[z]
   perform_collect_face(z, l1, h1, l2, h2, 3, ny[zone]#=-2=#, nz[zone]#=-2=#)
end

function collect_face(z, l1, h1, l2, h2, _::Val{4})
   zone = proc_zone_id[z]
   perform_collect_face(z, l1, h1, l2, h2, 4, ny[zone]#=-2=#, nz[zone]#=-2=#)
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

function go_cluster(clusters, niter, inorm, dt, ratio, x_zones, y_zones, gx_size, gy_size, gz_size, nxmax, nx_, ny_, nz_, proc_num_zones_all, x_size, y_size, zone_proc_id_, proc_zone_id_all, iz_west_, iz_east_, iz_south_, iz_north_, itimer, npb_verbose)
   global clusterid = Distributed.myid(role=:worker) - 2 

   global nx = nx_
   global ny = ny_
   global nz = nz_

   global zone_proc_id = zone_proc_id_
   global proc_zone_id = proc_zone_id_all[clusterid+1]
   global iz_west = iz_west_ 
   global iz_east = iz_east_ 
   global iz_south = iz_south_
   global iz_north = iz_north_

   global proc_num_zones = proc_num_zones_all[clusterid+1]

   global face_out = Array{Array{Array{Float64}}}(undef, proc_num_zones)
   for iz = 1:proc_num_zones
      zone = proc_zone_id[iz]
      face_out[iz] = Array{Array{Float64}}(undef, 4)
      face_out[iz][1] = Array{Float64}(undef, 5, ny[zone], gz_size) 
      face_out[iz][2] = Array{Float64}(undef, 5, ny[zone], gz_size) 
      face_out[iz][3] = Array{Float64}(undef, 5, nx[zone], gz_size) 
      face_out[iz][4] = Array{Float64}(undef, 5, nx[zone], gz_size) 
   end

   global face_in = Array{Array{Array{Float64}}}(undef, proc_num_zones)
   for iz = 1:proc_num_zones
      zone = proc_zone_id[iz]
      face_in[iz] = Array{Array{Float64}}(undef, 4)
      face_in[iz][1] = Array{Float64}(undef, 5, ny[zone], gz_size)
      face_in[iz][2] = Array{Float64}(undef, 5, ny[zone], gz_size)
      face_in[iz][3] = Array{Float64}(undef, 5, nx[zone], gz_size)
      face_in[iz][4] = Array{Float64}(undef, 5, nx[zone], gz_size)
   end

   global face_in_count = Array{Array{Int64}}(undef, proc_num_zones)
   global face_out_count = Array{Array{Int64}}(undef, proc_num_zones)
   global face_in_collect = Array{Array{Threads.Condition}}(undef, proc_num_zones)
   global face_in_receive = Array{Array{Threads.Condition}}(undef, proc_num_zones)
   for iz = 1:proc_num_zones
      zone = proc_zone_id[iz]
      face_in_count[iz] = Array{Int64}(undef, 4)
      face_out_count[iz] = zeros(Int64, 4)
      face_in_collect[iz] = Array{Threads.Condition}(undef, 4)
      face_in_receive[iz] = Array{Threads.Condition}(undef, 4)
      for f = 1:4
         face_in_collect[iz][f] = Threads.Condition()
         face_in_receive[iz][f] = Threads.Condition()
      end
      face_in_count[iz][1] = (ny[zone]#=-2=#)*(nz[zone]#=-2=#)*5
      face_in_count[iz][2] = (ny[zone]#=-2=#)*(nz[zone]#=-2=#)*5
      face_in_count[iz][3] = (nx[zone]#=-2=#)*(nz[zone]#=-2=#)*5 
      face_in_count[iz][4] = (nx[zone]#=-2=#)*(nz[zone]#=-2=#)*5
   end

   global lkout = ReentrantLock()
   global lkin = ReentrantLock()

   global lk1 = ReentrantLock()
   global lk2 = ReentrantLock()

   LU.go_node(clusterid, clusters, niter, inorm, dt, ratio, x_zones, y_zones, gx_size, gy_size, gz_size, nxmax, nx, ny, nz, proc_num_zones, proc_zone_id, iz_west, iz_east, iz_south, iz_north, itimer, npb_verbose) 
end

function go_node(clusterid, no_nodes, niter, inorm, dt, ratio, x_zones, y_zones, gx_size, gy_size, gz_size, nxmax, nx, ny, nz, proc_num_zones, proc_zone_id, iz_west, iz_east, iz_south, iz_north, itimer, npb_verbose)
   @everywhere workers() LU.perform($clusterid, $no_nodes, $niter, $inorm, $dt, $ratio, $x_zones, $y_zones, $gx_size, $gy_size, $gz_size, $nxmax, $nx, $ny, $nz, $proc_num_zones, $proc_zone_id, $iz_west, $iz_east, $iz_south, $iz_north, $itimer, $npb_verbose) 
end