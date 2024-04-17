
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

function exch_qbc(proc_num_zones, zone_proc_id, proc_zone_id, u, row, col, west, east, north, south, west2, east2, north2, south2, nx, ny, nz, timeron, ist, iend, jst, jend, iz_west, iz_east, iz_south, iz_north, comm_exch, buf_exch_w_out, buf_exch_e_out, buf_exch_n_out, buf_exch_s_out, buf_exch_w_in, buf_exch_e_in, buf_exch_n_in, buf_exch_s_in)

   r = Ref{Int64}(1)
   requests = Array{MPI.Request}(undef,8*proc_num_zones)
   for i = 1:8*proc_num_zones
      requests[i] = MPI.REQUEST_NULL
   end

   for iz = 1:proc_num_zones
      zone = proc_zone_id[iz]
      exchange_qbc_send(iz, zone_proc_id, u[iz], row[iz], col[iz], west[iz], east[iz], north[iz], south[iz], west2[iz], east2[iz], north2[iz], south2[iz], nx[iz], ny[iz], nz[iz], timeron, zone, ist[iz], iend[iz], jst[iz], jend[iz], iz_west[zone], iz_east[zone], iz_south[zone], iz_north[zone], requests, r, comm_exch, buf_exch_w_out[iz], buf_exch_e_out[iz], buf_exch_n_out[iz], buf_exch_s_out[iz])
   end


   for iz = 1:proc_num_zones
      zone = proc_zone_id[iz]
      west_zone = iz_west[zone]
      east_zone = iz_east[zone]
      north_zone = iz_south[zone]
      south_zone = iz_north[zone]
      zw = proc_zone_id_inv(west_zone, proc_zone_id) ; tag_west  = isnothing(zw) ? -1 : from_w + zw
      ze = proc_zone_id_inv(east_zone, proc_zone_id) ; tag_east  = isnothing(ze) ? -1 : from_e + ze
      zn = proc_zone_id_inv(north_zone, proc_zone_id); tag_north = isnothing(zn) ? -1 : from_n + zn
      zs = proc_zone_id_inv(south_zone, proc_zone_id); tag_south = isnothing(zs) ? -1 : from_s + zs

     #@info "$clusterid/$node: z=$iz --- RECV FROM NORTH src=$(north2[iz]) tag=$tag_north"
     #@info "$clusterid/$node: z=$iz --- RECV FROM SOUTH src=$(south2[iz]) tag=$tag_south"
     #@info "$clusterid/$node: z=$iz --- RECV FROM EAST src=$(east2[iz]) tag=$tag_east"
     #@info "$clusterid/$node: z=$iz --- RECV FROM WEST src=$(west2[iz]) tag=$tag_west"

      requests[r[]+0] = tag_west  > 0 && west[iz] == -1  && west_zone != -1  ? MPI.Irecv!(buf_exch_w_in[iz], comm_exch; source = west2[iz],  tag = tag_west) : MPI.REQUEST_NULL
      requests[r[]+1] = tag_east  > 0 && east[iz] == -1  && east_zone != -1  ? MPI.Irecv!(buf_exch_e_in[iz], comm_exch; source = east2[iz],  tag = tag_east) : MPI.REQUEST_NULL
      requests[r[]+2] = tag_north > 0 && north[iz] == -1 && north_zone != -1 ? MPI.Irecv!(buf_exch_n_in[iz], comm_exch; source = north2[iz], tag = tag_north) : MPI.REQUEST_NULL
      requests[r[]+3] = tag_south > 0 && south[iz] == -1 && south_zone != -1 ? MPI.Irecv!(buf_exch_s_in[iz], comm_exch; source = south2[iz], tag = tag_south) : MPI.REQUEST_NULL

      r[] = r[] + 4
   end

   if r[] > 1
       MPI.Waitall(requests)
   end
   

   for iz = 1:proc_num_zones
      zone = proc_zone_id[iz]
      exchange_qbc_recv(iz, zone_proc_id, proc_zone_id, u[iz], row[iz], col[iz], west[iz], east[iz], north[iz], south[iz], nx[iz], ny[iz], nz[iz], timeron, zone, ist[iz], iend[iz], jst[iz], jend[iz], iz_west[zone], iz_east[zone], iz_south[zone], iz_north[zone], comm_exch, buf_exch_w_in[iz], buf_exch_e_in[iz], buf_exch_n_in[iz], buf_exch_s_in[iz])
   end


end


function exchange_qbc_send(z, zone_proc_id, u, row, col, west, east, north, south, west2, east2, north2, south2, nx, ny, nz, timeron, zone, ist, iend, jst, jend, west_zone, east_zone, south_zone, north_zone, requests, r, comm_exch, buf_exch_w_out, buf_exch_e_out, buf_exch_n_out, buf_exch_s_out)

   #@info "$clusterid/$node --- COPY FACES SEND - BEGIN - z=$z zone=$zone"
                     
#---------------------------------------------------------------------
#     because the difference stencil for the diagonalized scheme is 
#     orthogonal, we do not have to perform the staged copying of faces, 
#     but can send all face information simultaneously to the neighboring 
#     cells in all directions          
#---------------------------------------------------------------------
      if (timeron) timer_start(t_rdis1) end

#---------------------------------------------------------------------
#     fill the buffer to be sent to northern neighbors (j_dir)
#---------------------------------------------------------------------
      if east == -1 && east_zone != -1
         u_face = view(u, 1:5, 1:nx, ny-1, 1:nz)
         zone_cluster = zone_proc_id[east_zone]
         if zone_cluster != clusterid
            ist = (row-1)*nx+1
            iend = row*nx
            kst = 1
            kend = nz
            #@info "$clusterid/$node: copy faces send EAST - BEGIN - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend " 
            remotecall(deposit_face, 1, z, ist, iend, kst, kend, u_face, Val(2); role=:worker)                                       
            #@info "$clusterid/$node: copy faces send EAST - END - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend"
         else
            tag_east = from_w + z 
            buf_exch_e_out .= reshape(u_face, 5*nx*nz)
           #@info "$clusterid/$node: z=$z --- SEND TO EAST dest=$east2 tag=$tag_east"
            requests[r[]] = MPI.Isend(buf_exch_e_out, comm_exch; dest = east2, tag = tag_east) 
            r[] += 1
         end
      end

#---------------------------------------------------------------------
#     fill the buffer to be sent to southern neighbors 
#---------------------------------------------------------------------
      if west == -1 && west_zone != -1
         u_face = view(u, 1:5, 1:nx, 2, 1:nz)
         zone_cluster = zone_proc_id[west_zone]
         if zone_cluster != clusterid
            ist = (row-1)*nx+1
            iend = row*nx
            kst = 1
            kend = nz
            #@info "$clusterid/$node: copy faces send WEST - BEGIN - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend " 
            remotecall(deposit_face, 1, z, ist, iend, kst, kend, u_face, Val(1); role=:worker)
            #@info "$clusterid/$node: copy faces send WEST - END - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend"
         else
            tag_west = from_e + z
            buf_exch_w_out .= reshape(u_face, 5*nx*nz)
           #@info "$clusterid/$node: z=$z --- SEND TO WEST dest=$west2 tag=$tag_west"
            requests[r[]] = MPI.Isend(buf_exch_w_out, comm_exch; dest = west2, tag = tag_west) 
            r[] += 1
         end
      end

#---------------------------------------------------------------------
#     fill the buffer to be sent to western neighbors 
#---------------------------------------------------------------------

      if north == -1 && north_zone != -1
         u_face = view(u, 1:5, 2, 1:ny, 1:nz)
         zone_cluster = zone_proc_id[north_zone]
         if zone_cluster != clusterid
            jst =  (col-1)*ny+1
            jend = col*ny
            kst = 1
            kend = nz
            #@info "$clusterid/$node: copy faces send NORTH - BEGIN - z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend " 
            remotecall(deposit_face, 1, z, jst, jend, kst, kend, u_face, Val(4); role=:worker)
            #@info "$clusterid/$node: copy faces send NORTH - END - z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend" 
         else
            tag_north = from_s + z
            buf_exch_n_out .= reshape(u_face, 5*ny*nz)
           #@info "$clusterid/$node: z=$z --- SEND TO NORTH dest=$north2 tag=$tag_north"
            requests[r[]] = MPI.Isend(buf_exch_n_out, comm_exch; dest = north2, tag = tag_north) 
            r[] += 1
         end
      end

#---------------------------------------------------------------------
#     fill the buffer to be sent to eastern neighbors (i-dir)
#---------------------------------------------------------------------
      if south == -1 && south_zone != -1
         u_face = view(u, 1:5, nx-1, 1:ny, 1:nz)
         zone_cluster = zone_proc_id[south_zone]
         if zone_cluster != clusterid
            jst =  (col-1)*ny+1
            jend = col*ny
            kst = 1
            kend = nz
            #@info "$clusterid/$node: copy faces send SOUTH - BEGIN - z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend " 
            remotecall(deposit_face, 1, z, jst, jend, kst, kend, u_face, Val(3); role=:worker)
            #@info "$clusterid/$node: copy faces send SOUTH - END - z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend" 
         else
            tag_south = from_n + z
            buf_exch_s_out .= reshape(u_face, 5*ny*nz)
           #@info "$clusterid/$node: z=$z --- SEND TO SOUTH dest=$south2 tag=$tag_south"
            requests[r[]] = MPI.Isend(buf_exch_s_out, comm_exch; dest = south2, tag = tag_south) 
            r[] += 1
         end
      end

       
      if (timeron) timer_stop(t_rdis1) end

   end

   function exchange_qbc_recv(z, zone_proc_id, proc_zone_id, u, row, col, west, east, north, south, nx, ny, nz, timeron, zone, ist, iend, jst, jend, west_zone, east_zone, south_zone, north_zone, comm_exch, buf_exch_w_in, buf_exch_e_in, buf_exch_n_in, buf_exch_s_in)

      if (timeron) timer_start(t_rdis2) end

#---------------------------------------------------------------------
#     unpack the data that has just been received;             
#---------------------------------------------------------------------

      if east == -1 && east_zone != -1
         zone_cluster = zone_proc_id[east_zone]
         if zone_cluster != clusterid
            ist = (row-1)*nx+1
            iend = row*nx
            kst = 1
            kend = nz
            #@info "$clusterid/$node: copy faces recv EAST - BEGIN - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend" 
            u_face = remotecall(collect_face, 1, z, ist, iend, kst, kend, Val(2); role=:worker)
            u_face_fetch = fetch(u_face; role=:worker)
            view(u, 1:5, 1:nx, ny, 1:nz) .= u_face_fetch
            #@info "$clusterid/$node: copy faces recv EAST - END - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend " 
         else            
            view(u, 1:5, 1:nx, ny, 1:nz) .= reshape(buf_exch_e_in, 5, nx, nz)
         end
      end

      if west == -1 && west_zone != -1
         zone_cluster = zone_proc_id[west_zone]
         if zone_cluster != clusterid
            ist = (row-1)*nx+1
            iend = row*nx
            kst = 1
            kend = nz
            #@info "$clusterid/$node: copy faces recv WEST - BEGIN - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend" 
            u_face = remotecall(collect_face, 1, z, ist, iend, kst, kend, Val(1); role=:worker)
            u_face_fetch = fetch(u_face; role=:worker)
            view(u, 1:5, 1:nx, 1, 1:nz) .= u_face_fetch
            #@info "$clusterid/$node: copy faces recv WEST - END - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend " 
         else
            view(u, 1:5, 1:nx, 1, 1:nz) .= reshape(buf_exch_w_in, 5, nx, nz)
         end
      end

      if north == -1 && north_zone != -1
         zone_cluster = zone_proc_id[north_zone]
         if zone_cluster != clusterid
            jst = (col-1)*ny+1
            jend = col*ny
            kst = 1
            kend = nz
            #@info "$clusterid/$node: copy faces recv NORTH - BEGIN - z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend" 
            u_face = remotecall(collect_face, 1, z, jst, jend, kst, kend, Val(4); role=:worker)
            u_face_fetch = fetch(u_face; role=:worker)
            view(u, 1:5, 1, 1:ny, 1:nz) .= u_face_fetch
            #@info "$clusterid/$node: copy faces recv NORTH - END - z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend " 
         else
            view(u, 1:5, 1, 1:ny, 1:nz) .= reshape(buf_exch_n_in, 5, ny, nz)
         end
      end

      if south == -1 && south_zone != -1
         zone_cluster = zone_proc_id[south_zone]
         if zone_cluster != clusterid
            jst = (col-1)*ny+1
            jend = col*ny
            kst = 1
            kend = nz
            #@info "$clusterid/$node: copy faces recv SOUTH - BEGIN - z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend" 
            u_face = remotecall(collect_face, 1, z, jst, jend, kst, kend, Val(3); role=:worker)
            u_face_fetch = fetch(u_face; role=:worker)
            view(u, 1:5, nx, 1:ny, 1:nz) .= u_face_fetch
            #@info "$clusterid/$node: copy faces recv SOUTH  - END- z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend " 
         else
            view(u, 1:5, nx, 1:ny, 1:nz) .= reshape(buf_exch_s_in, 5, ny, nz)
         end
      end

      if (timeron) timer_stop(t_rdis2) end

      return nothing
end
