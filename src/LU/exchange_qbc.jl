
#---------------------------------------------------------------------
# This function copies the face values of a variable defined on a set 
# of cells to the overlap locations of the adjacent sets of cells. 
# Because a set of cells interfaces in each direction with exactly one 
# other set, we only need to fill six different buffers. We could try to 
# overlap communication with computation, by computing
# some internal values while communicating boundary values, but this
# adds so much overhead that it's not clearly useful. 
#---------------------------------------------------------------------


#---------------------------------------------------------------------
#     exit immediately if there are no faces to be copied           
#---------------------------------------------------------------------

function exch_qbc(proc_num_zones, u, row, col, west, east, north, south, nx, ny, nz, timeron, proc_zone_id, ist, iend, jst, jend, iz_west, iz_east, iz_south, iz_north)

   for iz = 1:proc_num_zones
      zone = proc_zone_id[iz]
      copy_faces_send(iz, u[iz], row[iz], col[iz], west[iz], east[iz], north[iz], south[iz], nx[iz], ny[iz], nz[iz], timeron, zone, ist[iz], iend[iz], jst[iz], jend[iz], iz_west[zone], iz_east[zone], iz_south[zone], iz_north[zone])
   end

   for iz = 1:proc_num_zones
      zone = proc_zone_id[iz]
      copy_faces_recv(iz, u[iz], row[iz], col[iz], west[iz], east[iz], north[iz], south[iz], nx[iz], ny[iz], nz[iz], timeron, zone, ist[iz], iend[iz], jst[iz], jend[iz], iz_west[zone], iz_east[zone], iz_south[zone], iz_north[zone])
   end

end


function copy_faces_send(z, u, row, col, west, east, north, south, nx, ny, nz, timeron, zone, ist, iend, jst, jend, west_zone, east_zone, south_zone, north_zone)

   #@info "$clusterid/$node --- COPY FACES SEND - BEGIN - z=$z zone=$zone"
                     
#---------------------------------------------------------------------
#     because the difference stencil for the diagonalized scheme is 
#     orthogonal, we do not have to perform the staged copying of faces, 
#     but can send all face information simultaneously to the neighboring 
#     cells in all directions          
#---------------------------------------------------------------------
#      if (timeron) timer_start(t_rdis1) end


#---------------------------------------------------------------------
#     fill the buffer to be sent to northern neighbors (j_dir)
#---------------------------------------------------------------------
      if east == -1 && east_zone != -1
         #a = north == -1 ? 1 : 0
         #b = south == -1 ? 1 : 0
         nnx = nx #- a - b
         u_face = view(u, 1:5, 1:nnx, ny-1, 1:nz#=-2=#)
         ist = (row-1)*nnx+1
         iend = row*nnx
         kst = 1
         kend = nz#=-2=#
         #@info "$clusterid/$node: copy faces send EAST - BEGIN - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend -- u_face = $u_face" 
         remotecall(deposit_face, 1, z, ist, iend, kst, kend, u_face, Val(2); role=:worker)                                       
         #@info "$clusterid/$node: copy faces send EAST - END - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend"
      end

#---------------------------------------------------------------------
#     fill the buffer to be sent to southern neighbors 
#---------------------------------------------------------------------
      if west == -1 && west_zone != -1
         #a = north == -1 ? 1 : 0
         #b = south == -1 ? 1 : 0
         nnx = nx #- a - b
         u_face = view(u, 1:5, 1:nnx, 2, 1:nz#=-2=#)
         ist = (row-1)*nnx+1
         iend = row*nnx
         kst = 1
         kend = nz#=-2=#
         #@info "$clusterid/$node: copy faces send WEST - BEGIN - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend -- u_face = $u_face" 
         remotecall(deposit_face, 1, z, ist, iend, kst, kend, u_face, Val(1); role=:worker)
         #@info "$clusterid/$node: copy faces send WEST - END - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend"
      end

#---------------------------------------------------------------------
#     fill the buffer to be sent to western neighbors 
#---------------------------------------------------------------------

      if north == -1 && north_zone != -1
         #a = west == -1 ? 1 : 0
         #b = east == -1 ? 1 : 0
         nny = ny #- a - b
         u_face = view(u, 1:5, 2, 1:nny, 1:nz#=-2=#)
         jst =  (col-1)*nny+1
         jend = col*nny
         kst = 1
         kend = nz#=-2=#
         #@info "$clusterid/$node: copy faces send NORTH - BEGIN - z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend -- u_face = $u_face" 
         remotecall(deposit_face, 1, z, jst, jend, kst, kend, u_face, Val(4); role=:worker)
         #@info "$clusterid/$node: copy faces send NORTH - END - z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend" 
      end

#---------------------------------------------------------------------
#     fill the buffer to be sent to eastern neighbors (i-dir)
#---------------------------------------------------------------------
      if south == -1 && south_zone != -1
         #a = west == -1 ? 1 : 0
         #b = east == -1 ? 1 : 0
         nny = ny #- a - b
         u_face = view(u, 1:5, nx-1, 1:nny, 1:nz#=-2=#)
         jst =  (col-1)*nny+1
         jend = col*nny
         kst = 1
         kend = nz#=-2=#
         #@info "$clusterid/$node: copy faces send SOUTH - BEGIN - z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend -- u_face = $u_face" 
         remotecall(deposit_face, 1, z, jst, jend, kst, kend, u_face, Val(3); role=:worker)
         #@info "$clusterid/$node: copy faces send SOUTH - END - z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend" 
      end


      #@info "$clusterid/$node --- COPY FACES SEND - END - z=$z zone=$zone"

 #     if (timeron) timer_stop(t_rdis1) end

   end

   function copy_faces_recv(z, u, row, col, west, east, north, south, nx, ny, nz, timeron, zone, ist, iend, jst, jend, west_zone, east_zone, south_zone, north_zone)

      #@info "$clusterid/$node --- COPY FACES RECV - BEGIN - z=$z zone=$zone "

#      if (timeron) timer_start(t_rdis2) end

#---------------------------------------------------------------------
#     unpack the data that has just been received;             
#---------------------------------------------------------------------

      if east == -1 && east_zone != -1
         #a = north == -1 ? 1 : 0
         #b = south == -1 ? 1 : 0
         nnx = nx #- a - b
         ist = (row-1)*nnx+1
         iend = row*nnx
         kst = 1
         kend = nz#=-2=#
         #@info "$clusterid/$node: copy faces recv EAST - BEGIN - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend" 
         u_face = remotecall(collect_face, 1, z, ist, iend, kst, kend, Val(2); role=:worker)
         u_face_fetch = fetch(u_face; role=:worker)
         view(u, 1:5, 1:nnx, ny, 1:nz#=-2=#) .= u_face_fetch
         #@info "$clusterid/$node: copy faces recv EAST - END - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend -- u_face = $u_face_fetch" 
      end

      if west == -1 && west_zone != -1
         #a = north == -1 ? 1 : 0
         #b = south == -1 ? 1 : 0
         nnx = nx #- a - b
         ist = (row-1)*nnx+1
         iend = row*nnx
         kst = 1
         kend = nz#=-2=#
         #@info "$clusterid/$node: copy faces recv WEST - BEGIN - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend" 
         u_face = remotecall(collect_face, 1, z, ist, iend, kst, kend, Val(1); role=:worker)
         u_face_fetch = fetch(u_face; role=:worker)
         view(u, 1:5, 1:nnx, 1, 1:nz#=-2=#) .= u_face_fetch
         #@info "$clusterid/$node: copy faces recv WEST - END - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend -- u_face = $u_face_fetch" 
      end

      if north == -1 && north_zone != -1
         #a = west == -1 ? 1 : 0
         #b = east == -1 ? 1 : 0
         nny = ny #- a - b
         jst = (col-1)*nny+1
         jend = col*nny
         kst = 1
         kend = nz#=-2=#
         #@info "$clusterid/$node: copy faces recv NORTH - BEGIN - z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend" 
         u_face = remotecall(collect_face, 1, z, jst, jend, kst, kend, Val(4); role=:worker)
         u_face_fetch = fetch(u_face; role=:worker)
         view(u, 1:5, 1, 1:nny, 1:nz#=-2=#) .= u_face_fetch
         #@info "$clusterid/$node: copy faces recv NORTH - END - z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend -- u_face = $u_face_fetch" 
      end

      if south == -1 && south_zone != -1
         #a = west == -1 ? 1 : 0
         #b = east == -1 ? 1 : 0
         nny = ny #- a - b
         jst = (col-1)*nny+1
         jend = col*nny
         kst = 1
         kend = nz#=-2=#
         #@info "$clusterid/$node: copy faces recv SOUTH - BEGIN - z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend" 
         u_face = remotecall(collect_face, 1, z, jst, jend, kst, kend, Val(3); role=:worker)
         u_face_fetch = fetch(u_face; role=:worker)
         view(u, 1:5, nx, 1:nny, 1:nz#=-2=#) .= u_face_fetch
         #@info "$clusterid/$node: copy faces recv SOUTH  - END- z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend -- u_face = $u_face_fetch" 
      end

      #@info "$clusterid/$node --- COPY FACES RECV - END - z=$z zone=$zone"

#      if (timeron) timer_stop(t_rdis2) end

      return nothing
end
