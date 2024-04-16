
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

function exch_qbc(proc_num_zones, zone_proc_id, proc_zone_id, u, row, col, west, east, north, south, nx, ny, nz, timeron, ist, iend, jst, jend, iz_west, iz_east, iz_south, iz_north, comm_exch, buf_exch_w_out, buf_exch_e_out, buf_exch_n_out, buf_exch_s_out, buf_exch_w_in, buf_exch_e_in, buf_exch_n_in, buf_exch_s_in)

   req_count = Ref{Int64}(1)
   requests = Array{MPI.Request}(undef,8*proc_num_zones)

   for iz = 1:proc_num_zones
      zone = proc_zone_id[iz]
      exchange_qbc_send(iz, zone_proc_id, u[iz], row[iz], col[iz], west[iz], east[iz], north[iz], south[iz], nx[iz], ny[iz], nz[iz], timeron, zone, ist[iz], iend[iz], jst[iz], jend[iz], iz_west[zone], iz_east[zone], iz_south[zone], iz_north[zone], requests, req_count, comm_exch, buf_exch_w_out[iz], buf_exch_e_out[iz], buf_exch_n_out[iz], buf_exch_s_out[iz])
   end

   for iz = 1:proc_num_zones
      zone = proc_zone_id[iz]
      exchange_qbc_recv(iz, zone_proc_id, proc_zone_id, u[iz], row[iz], col[iz], west[iz], east[iz], north[iz], south[iz], nx[iz], ny[iz], nz[iz], timeron, zone, ist[iz], iend[iz], jst[iz], jend[iz], iz_west[zone], iz_east[zone], iz_south[zone], iz_north[zone], requests, req_count, comm_exch, buf_exch_w_in[iz], buf_exch_e_in[iz], buf_exch_n_in[iz], buf_exch_s_in[iz])
   end

   if req_count[] > 1
      @info "$clusterid/$node: WAIT ALL BEGIN ! $requests"
      MPI.Waitall(requests)
      @info "$clusterid/$node: WAIT ALL END !"
   end

end


function exchange_qbc_send(z, zone_proc_id, u, row, col, west, east, north, south, nx, ny, nz, timeron, zone, ist, iend, jst, jend, west_zone, east_zone, south_zone, north_zone, requests, r, comm_exch, buf_exch_w_out, buf_exch_e_out, buf_exch_n_out, buf_exch_s_out)

   #@info "$clusterid/$node --- COPY FACES SEND - BEGIN - z=$z zone=$zone"
                     
#---------------------------------------------------------------------
#     because the difference stencil for the diagonalized scheme is 
#     orthogonal, we do not have to perform the staged copying of faces, 
#     but can send all face information simultaneously to the neighboring 
#     cells in all directions          
#---------------------------------------------------------------------
      if (timeron) timer_start(t_rdis1) end

      @info "$clusterid/$node: z=$z EXCHANGE_QBC SEND BEGIN"

#---------------------------------------------------------------------
#     fill the buffer to be sent to northern neighbors (j_dir)
#---------------------------------------------------------------------
      requests[r[]] = MPI.REQUEST_NULL
      if east == -1 && east_zone != -1
         u_face = view(u, 1:5, 1:nx, ny-1, 1:nz)
         zone_cluster = zone_proc_id[east_zone]
         if zone_cluster != clusterid
            @info "$clusterid/$node: z=$z SEND OUT 1 BEGIN"
            ist = (row-1)*nx+1
            iend = row*nx
            kst = 1
            kend = nz
            #@info "$clusterid/$node: copy faces send EAST - BEGIN - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend " 
            remotecall(deposit_face, 1, z, ist, iend, kst, kend, u_face, Val(2); role=:worker)                                       
            #@info "$clusterid/$node: copy faces send EAST - END - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend"
            @info "$clusterid/$node: z=$z SEND OUT 1 END"
         else
            @info "$clusterid/$node: z=$z SEND IN 1 BEGIN"
             #@info "$clusterid/$node: tag = $tag_west --- dest=$(predecessor[1]) -- z=$z SEND TO WEST -- b_size[2]=$(b_size[2]) -- ss[2]=$(ss[2])"
            tag_east = from_w + z 
            buf_exch_e_out .= reshape(u_face, 5*nx*nz)
            requests[r[]] = MPI.Isend(buf_exch_e_out, comm_exch; dest = east, tag = tag_east) 
            @info "$clusterid/$node: z=$z SEND IN 1 END"
         end
      end
      r[] += 1

      @info "$clusterid/$node: z=$z EXCHANGE_QBC SEND 1"

#---------------------------------------------------------------------
#     fill the buffer to be sent to southern neighbors 
#---------------------------------------------------------------------
      requests[r[]] = MPI.REQUEST_NULL
      if west == -1 && west_zone != -1
         u_face = view(u, 1:5, 1:nx, 2, 1:nz)
         zone_cluster = zone_proc_id[west_zone]
         if zone_cluster != clusterid
            @info "$clusterid/$node: z=$z SEND OUT 2 BEGIN"
            ist = (row-1)*nx+1
            iend = row*nx
            kst = 1
            kend = nz
            #@info "$clusterid/$node: copy faces send WEST - BEGIN - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend " 
            remotecall(deposit_face, 1, z, ist, iend, kst, kend, u_face, Val(1); role=:worker)
            #@info "$clusterid/$node: copy faces send WEST - END - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend"
            @info "$clusterid/$node: z=$z SEND OUT 2 END"
         else
            @info "$clusterid/$node: z=$z SEND IN 2 BEGIN"
            tag_west = from_e + z
            #@info "$clusterid/$node: tag = $tag_east --- dest=$(successor[1]) -- z=$z SEND TO EAST -- b_size[1]=$(b_size[1]) -- ss[1]=$(ss[1])"
            buf_exch_w_out .= reshape(u_face, 5*nx*nz)
            requests[r[]] = MPI.Isend(buf_exch_w_out, comm_exch; dest = west, tag = tag_west) 
            @info "$clusterid/$node: z=$z SEND IN 2 END"
         end
      end
      r[] += 1

      @info "$clusterid/$node: z=$z EXCHANGE_QBC SEND 2"

#---------------------------------------------------------------------
#     fill the buffer to be sent to western neighbors 
#---------------------------------------------------------------------

      requests[r[]] = MPI.REQUEST_NULL
      if north == -1 && north_zone != -1
         u_face = view(u, 1:5, 2, 1:ny, 1:nz)
         zone_cluster = zone_proc_id[north_zone]
         if zone_cluster != clusterid
            @info "$clusterid/$node: z=$z SEND OUT 3 BEGIN"
            jst =  (col-1)*ny+1
            jend = col*ny
            kst = 1
            kend = nz
            #@info "$clusterid/$node: copy faces send NORTH - BEGIN - z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend " 
            remotecall(deposit_face, 1, z, jst, jend, kst, kend, u_face, Val(4); role=:worker)
            #@info "$clusterid/$node: copy faces send NORTH - END - z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend" 
            @info "$clusterid/$node: z=$z SEND OUT 3 END"
         else
            @info "$clusterid/$node: z=$z SEND IN 3 BEGIN"
             #@info "$clusterid/$node: tag = $tag_south --- dest=$(successor[2]) -- z=$z SEND TO SOUTH -- b_size[3]=$(b_size[3]) -- ss[3]=$(ss[3])"
            tag_south = from_s + z
            buf_exch_n_out .= reshape(u_face, 5*ny*nz)
            requests[r[]] = MPI.Isend(buf_exch_n_out, comm_exch; dest = south, tag = tag_south) 
            @info "$clusterid/$node: z=$z SEND IN 3 END"
         end
      end
      r[] += 1

      @info "$clusterid/$node: z=$z EXCHANGE_QBC SEND 3"

#---------------------------------------------------------------------
#     fill the buffer to be sent to eastern neighbors (i-dir)
#---------------------------------------------------------------------
      requests[r[]] = MPI.REQUEST_NULL
      if south == -1 && south_zone != -1
         @info "$clusterid/$node: z=$z SEND OUT 4 BEGIN"
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
             @info "$clusterid/$node: z=$z SEND OUT 4 END"
         else
            @info "$clusterid/$node: z=$z SEND IN 4 BEGIN"
            #@info "$clusterid/$node: tag = $tag_north --- dest=$(predecessor[2]) -- z=$z SEND TO NORTH -- b_size[4]=$(b_size[4]) -- ss[4]=$(ss[4])"
            tag_north = from_n + z
            buf_exch_s_out .= reshape(u_face, 5*ny*nz)
            requests[r[]] = MPI.Isend(buf_exch_s_out, comm_exch; dest = north, tag = tag_north) 
            @info "$clusterid/$node: z=$z SEND IN 4 END"
         end
      end
      r[] += 1

       
      @info "$clusterid/$node: z=$z EXCHANGE_QBC SEND END"

      if (timeron) timer_stop(t_rdis1) end

   end

   function exchange_qbc_recv(z, zone_proc_id, proc_zone_id, u, row, col, west, east, north, south, nx, ny, nz, timeron, zone, ist, iend, jst, jend, west_zone, east_zone, south_zone, north_zone, requests, r, comm_exch, buf_exch_w_in, buf_exch_e_in, buf_exch_n_in, buf_exch_s_in)

      @info "$clusterid/$node: z=$z EXCHANGE_QBC RECV BEGIN"

      if (timeron) timer_start(t_rdis2) end

#---------------------------------------------------------------------
#     unpack the data that has just been received;             
#---------------------------------------------------------------------

      requests[r[]] = MPI.REQUEST_NULL
      if east == -1 && east_zone != -1
         zone_cluster = zone_proc_id[east_zone]
         if zone_cluster != clusterid
            @info "$clusterid/$node: z=$z RECV OUT 1 BEGIN"
            ist = (row-1)*nx+1
            iend = row*nx
            kst = 1
            kend = nz
            #@info "$clusterid/$node: copy faces recv EAST - BEGIN - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend" 
            u_face = remotecall(collect_face, 1, z, ist, iend, kst, kend, Val(2); role=:worker)
            u_face_fetch = fetch(u_face; role=:worker)
            view(u, 1:5, 1:nx, ny, 1:nz) .= u_face_fetch
            #@info "$clusterid/$node: copy faces recv EAST - END - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend " 
            @info "$clusterid/$node: z=$z RECV OUT 1 END"
         else
            @info "$clusterid/$node: z=$z RECV IN 1 BEGIN"
            ze = proc_zone_id_inv(east_zone, proc_zone_id); tag_east = from_e + ze
            requests[r[]] = MPI.Irecv!(buf_exch_e_in, comm_exch; source = east, tag = tag_east)
            view(u, 1:5, 1:nx, ny, 1:nz) .= reshape(buf_exch_e_in, 5, nx, nz)
            @info "$clusterid/$node: z=$z RECV IN 1 END"
         end
      end
      r[] += 1

      @info "$clusterid/$node: z=$z EXCHANGE_QBC RECV 1"

      requests[r[]] = MPI.REQUEST_NULL
      if west == -1 && west_zone != -1
         zone_cluster = zone_proc_id[west_zone]
         if zone_cluster != clusterid
            @info "$clusterid/$node: z=$z RECV OUT 2 BEGIN"
            ist = (row-1)*nx+1
            iend = row*nx
            kst = 1
            kend = nz
            #@info "$clusterid/$node: copy faces recv WEST - BEGIN - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend" 
            u_face = remotecall(collect_face, 1, z, ist, iend, kst, kend, Val(1); role=:worker)
            u_face_fetch = fetch(u_face; role=:worker)
            view(u, 1:5, 1:nx, 1, 1:nz) .= u_face_fetch
            #@info "$clusterid/$node: copy faces recv WEST - END - z=$z zone=$zone  ist=$ist  iend=$iend  kst=$kst  kend=$kend " 
             @info "$clusterid/$node: z=$z RECV OUT 2 END"
         else
            @info "$clusterid/$node: z=$z RECV IN 2 BEGIN 0 --- $(proc_zone_id_inv(west_zone, proc_zone_id)) zone_cluster=$zone_cluster ---proc_zone_id=$proc_zone_id zone_proc_id=$zone_proc_id --- west_zone=$west_zone zone=$zone"
            zw = proc_zone_id_inv(west_zone, proc_zone_id); tag_west = from_w + zw
            @info "$clusterid/$node: z=$z RECV IN 2 BEGIN 1"
            requests[r[]] = MPI.Irecv!(buf_exch_w_in, comm_exch; source = west, tag = tag_west) 
            @info "$clusterid/$node: z=$z RECV IN 2 BEGIN 2"
            view(u, 1:5, 1:nx, 1, 1:nz) .= reshape(buf_exch_e_in, 5, nx, nz)
            @info "$clusterid/$node: z=$z RECV IN 2 END"
         end
      end
      r[] += 1

      @info "$clusterid/$node: z=$z EXCHANGE_QBC RECV 2"

      requests[r[]] = MPI.REQUEST_NULL
      if north == -1 && north_zone != -1
         zone_cluster = zone_proc_id[north_zone]
         if zone_cluster != clusterid
            @info "$clusterid/$node: z=$z RECV OUT 3 BEGIN"
            jst = (col-1)*ny+1
            jend = col*ny
            kst = 1
            kend = nz
            #@info "$clusterid/$node: copy faces recv NORTH - BEGIN - z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend" 
            u_face = remotecall(collect_face, 1, z, jst, jend, kst, kend, Val(4); role=:worker)
            u_face_fetch = fetch(u_face; role=:worker)
            view(u, 1:5, 1, 1:ny, 1:nz) .= u_face_fetch
            #@info "$clusterid/$node: copy faces recv NORTH - END - z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend " 
            @info "$clusterid/$node: z=$z RECV OUT 3 END"
         else
            @info "$clusterid/$node: z=$z RECV IN 3 BEGIN"
            zn = proc_zone_id_inv(north_zone, proc_zone_id); tag_north = from_n + zn
            requests[r[]] = MPI.Irecv!(buf_exch_n_in, comm_exch; source = south, tag = tag_north) 
            view(u, 1:5, 1, 1:ny, 1:nz) .= reshape(buf_exch_n_in, 5, ny, nz)
            @info "$clusterid/$node: z=$z RECV IN 3 END"
         end
      end
      r[] += 1

      @info "$clusterid/$node: z=$z EXCHANGE_QBC RECV 3"

      requests[r[]] = MPI.REQUEST_NULL
      if south == -1 && south_zone != -1
         zone_cluster = zone_proc_id[south_zone]
         if zone_cluster != clusterid
            @info "$clusterid/$node: z=$z RECV OUT 4 BEGIN"
            jst = (col-1)*ny+1
            jend = col*ny
            kst = 1
            kend = nz
            #@info "$clusterid/$node: copy faces recv SOUTH - BEGIN - z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend" 
            u_face = remotecall(collect_face, 1, z, jst, jend, kst, kend, Val(3); role=:worker)
            u_face_fetch = fetch(u_face; role=:worker)
            view(u, 1:5, nx, 1:ny, 1:nz) .= u_face_fetch
            #@info "$clusterid/$node: copy faces recv SOUTH  - END- z=$z zone=$zone  jst=$jst  jend=$jend  kst=$kst  kend=$kend " 
            @info "$clusterid/$node: z=$z RECV OUT 4 END"
         else
            @info "$clusterid/$node: z=$z RECV IN 4 BEGIN"
            zs = proc_zone_id_inv(south_zone, proc_zone_id); tag_south = from_s + zs
            requests[r[]] = MPI.Irecv!(buf_exch_s_in, comm_exch; source = north, tag = tag_south)
            view(u, 1:5, 1, 1:ny, 1:nz) .= reshape(buf_exch_s_in, 5, ny, nz)
            @info "$clusterid/$node: z=$z RECV IN 4 END"
         end
      end
      r[] += 1

      @info "$clusterid/$node: z=$z EXCHANGE_QBC RECV END"

      if (timeron) timer_stop(t_rdis2) end

      return nothing
end
