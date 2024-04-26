function exch_qbc(_::Val{ncells},
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
                    timeron,
                    ) where {ncells}  

    for iz = 1:proc_num_zones
        exch_qbc_send(Val(ncells),
                        iz,
                        zone_proc_id,
                        proc_zone_id,
                        iz_west,
                        iz_east,
                        iz_north,
                        iz_south,
                        cell_coord[iz],
                        cell_size[iz],
                        cell_start[iz],
                        cell_end[iz],
                        cell_low[iz],
                        cell_high[iz],
                        successor[iz],
                        predecessor[iz],
                        u[iz],
                        in_buffer[iz],
                        out_buffer[iz],
                        requests[iz],
                        ss[iz],
                        sr[iz],
                        b_size[iz],
                        comm_exch[iz],
                        timeron,)
    end
 
    for iz = 1:proc_num_zones
        exch_qbc_recv(Val(ncells),
                        iz,
                        zone_proc_id,
                        proc_zone_id,
                        iz_west,
                        iz_east,
                        iz_north,
                        iz_south,
                        cell_coord[iz],
                        cell_size[iz],
                        cell_start[iz],
                        cell_end[iz],
                        cell_low[iz],
                        cell_high[iz],
                        successor[iz],
                        predecessor[iz],
                        u[iz], 
                        in_buffer[iz],
                        out_buffer[iz],
                        requests[iz],
                        ss[iz],
                        sr[iz],
                        b_size[iz],
                        comm_exch[iz],
                        timeron,)
    end

    for iz = 1:proc_num_zones
    end

 end
 

#---------------------------------------------------------------------
# this function copies the face values of a variable defined on a set 
# of cells to the overlap locations of the adjacent sets of cells. 
# Because a set of cells interfaces in each direction with exactly one 
# other set, we only need to fill six different buffers. We could try to 
# overlap communication with computation, by computing
# some internal values while communicating boundary values, but this
# adds so much overhead that it's not clearly useful. 
#---------------------------------------------------------------------

const TO_WEST = 1
const TO_EAST = 2
const TO_SOUTH = 3
const TO_NORTH = 4

const FROM_WEST = 1
const FROM_EAST = 2
const FROM_SOUTH = 3
const FROM_NORTH = 4

function exch_qbc_send(_::Val{ncells},
                        z,
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
                        timeron,
                        ) where {ncells}  
                           
#---------------------------------------------------------------------
# because the difference stencil for the diagonalized scheme is 
# orthogonal, we do not have to perform the staged copying of faces, 
# but can send all face information simultaneously to the neighboring 
# cells in all directions          
#---------------------------------------------------------------------
       if (timeron) timer_start(t_rdis1) end
 
       p0 = 0
       p1 = 0
       p2 = 0
       p3 = 0
 
       for c = 1:ncells

#---------------------------------------------------------------------
#            fill the buffer to be sent to western neighbors 
#---------------------------------------------------------------------
             if cell_coord[1, c] == ncells
                # outer face (inter-zone)
               zone_cluster = zone_proc_id[iz_east[proc_zone_id[z]]]
               if zone_cluster != clusterid
                  u_face = view(u, cell_size[1, c]-2, cell_start[2,c]:cell_size[2, c]-cell_end[2,c]-1, 
                                                      cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, 1:5, c)
                  remotecall(deposit_face, 1, z, cell_low[2, c] + cell_start[2,c] + 1, cell_high[2, c] - cell_end[2,c] + 1, 
                                                 cell_low[3, c] + cell_start[3,c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                                 u_face, Val(TO_EAST); role=:worker)
               else
                  #@info "$clusterid/$node: zone=$(proc_zone_id[z]) other_zone=$(iz_east[proc_zone_id[z]]) zone_cluster=$zone_cluster TO_EAST"
                  for m = 1:5
                      for k = cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1
                          for j = cell_start[2,c]:cell_size[2, c]-cell_end[2,c]-1
                              i = cell_size[1, c]-2
                            @info "$clusterid/$node: z=$z send-u[$i,$j,$k,$m,$c]=$(u[i,j,k,m,c]) - TO EAST p0=$p0 ss[1]=$(ss[1])"
                              out_buffer[ss[1] + p0] = u[i, j, k, m, c]
                              p0 = p0 + 1
                          end
                      end
                  end
               end
             end

#---------------------------------------------------------------------
#            fill the buffer to be sent to eastern neighbors (i-dir)
#---------------------------------------------------------------------
            if cell_coord[1, c] == 1
               # outer face (inter-zone)
               zone_cluster = zone_proc_id[iz_west[proc_zone_id[z]]]
               if zone_cluster != clusterid
                  # adjacent zone is in another cluster (send through Distribute.jl)
                  u_face = view(u, 1, cell_start[2,c]:cell_size[2, c]-cell_end[2,c]-1, 
                                      cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, 1:5, c)
                  # TO WEST
                  remotecall(deposit_face, 1, z, cell_low[2, c] + cell_start[2,c] + 1, cell_high[2, c] - cell_end[2,c] + 1, 
                                                 cell_low[3, c] + cell_start[3,c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                                 u_face, Val(TO_WEST); role=:worker)
               else
                  #@info "$clusterid/$node: zone=$(proc_zone_id[z]) other_zone=$(iz_west[proc_zone_id[z]]) zone_cluster=$zone_cluster TO_WEST"
                  # adjacent zone is the same cluster (send through MPI.jl)
                  for m = 1:5
                      for k = cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1
                          for j = cell_start[2,c]:cell_size[2, c]-cell_end[2,c]-1
                              i = 1
                            @info "$clusterid/$node: z=$z send-u[$i,$j,$k,$m,$c]=$(u[i,j,k,m,c]) - TO WEST p1=$p1 ss[2]=$(ss[2])"
                              out_buffer[ss[2] + p1] = u[i, j, k, m, c]
                              p1 = p1 + 1
                          end
                        end
                  end
               end
            end

#---------------------------------------------------------------------
#            fill the buffer to be sent to northern neighbors (j_dir)
#---------------------------------------------------------------------
            if cell_coord[2, c] == ncells
                # outer face (inter-zone)
               zone_cluster = zone_proc_id[iz_south[proc_zone_id[z]]]
               if zone_cluster != clusterid
                  u_face = view(u, cell_start[1,c]:cell_size[1, c]-cell_end[1,c]-1, cell_size[2, c]-2, 
                                   cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, 1:5, c)
                  remotecall(deposit_face, 1, z, cell_low[1, c] + cell_start[1,c] + 1, cell_high[1, c] - cell_end[1,c] + 1, 
                                                 cell_low[3, c] + cell_start[3,c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                                 u_face, Val(TO_SOUTH); role=:worker)
               else
                  #@info "$clusterid/$node: zone=$(proc_zone_id[z]) other_zone=$(iz_south[proc_zone_id[z]]) zone_cluster=$zone_cluster TO_SOUTH"
                  for m = 1:5
                      for k = cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1
                          j = cell_size[2, c]-2
                          for i = cell_start[1,c]:cell_size[1, c]-cell_end[1,c]-1
                            @info "$clusterid/$node: z=$z send-u[$i,$j,$k,$m,$c]=$(u[i,j,k,m,c]) - TO SOUTH p2=$p2 ss[3]=$(ss[3])"
                              out_buffer[ss[3] + p2] = u[i, j, k, m, c]
                              p2 = p2 + 1
                        end
                      end
                  end
               end
            end

#---------------------------------------------------------------------
#            fill the buffer to be sent to southern neighbors 
#---------------------------------------------------------------------
             if cell_coord[2, c] == 1
                u_face = view(u, cell_start[1,c]:cell_size[1, c]-cell_end[1,c]-1, 1, 
                                 cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, 1:5, c)
                zone_cluster = zone_proc_id[iz_north[proc_zone_id[z]]]
                if zone_cluster != clusterid
                   remotecall(deposit_face, 1, z, cell_low[1, c] + cell_start[1,c] + 1, cell_high[1, c] - cell_end[1,c] + 1, 
                                                  cell_low[3, c] + cell_start[3,c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                                  u_face, Val(TO_NORTH); role=:worker)
                else
                   #@info "$clusterid/$node: zone=$(proc_zone_id[z]) other_zone=$(iz_north[proc_zone_id[z]]) zone_cluster=$zone_cluster TO_NORTH"
                   for m = 1:5
                       for k = cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1
                           j = 1    
                            for i = cell_start[1,c]:cell_size[1, c]-cell_end[1,c]-1
                              @info "$clusterid/$node: z=$z send-u[$i,$j,$k,$m,$c]=$(u[i,j,k,m,c]) - TO NORTH p3=$p3 ss[4]=$(ss[4])"
                                out_buffer[ss[4] + p3] = u[i, j, k, m, c]
                                p3 = p3 + 1
                           end
                       end
                   end
                end
            end


#---------------------------------------------------------------------
#       cell loop
#---------------------------------------------------------------------
       end

     @info "$clusterid/$node: z=$z SEND TO EAST --- dest=$(successor[1]) -- tag = $WEST -- b_size[1]=$(b_size[1])"
     @info "$clusterid/$node: z=$z SEND TO WEST --- dest=$(predecessor[1]) -- tag = $EAST -- b_size[2]=$(b_size[2])"
     @info "$clusterid/$node: z=$z SEND TO SOUTH --- dest=$(successor[2]) -- tag = $SOUTH -- b_size[3]=$(b_size[3])"
     @info "$clusterid/$node: z=$z SEND TO NORTH --- dest=$(predecessor[2]) -- tag = $NORTH -- b_size[4]=$(b_size[4])"

       if z==5 && clusterid==0 && node==2
           for i = eachindex(out_buffer)
          @info "$clusterid/$node: z=$z out_buffer[$i] = $(out_buffer[i])"
           end
       end

      requests[5] = MPI.Isend(view(out_buffer, ss[1]:ss[1]+b_size[1]-1), comm_exch; dest = successor[1], tag = WEST)
      requests[6] = MPI.Isend(view(out_buffer, ss[2]:ss[2]+b_size[2]-1), comm_exch; dest = predecessor[1], tag = EAST)
      requests[7] = MPI.Isend(view(out_buffer, ss[3]:ss[3]+b_size[3]-1), comm_exch; dest = successor[2], tag = SOUTH)
      requests[8] = MPI.Isend(view(out_buffer, ss[4]:ss[4]+b_size[4]-1), comm_exch; dest = predecessor[2], tag = NORTH)

       if (timeron) timer_stop(t_rdis1) end

end


function exch_qbc_recv(_::Val{ncells},
                        z,
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
                        timeron,
                     ) where {ncells}  

#---------------------------------------------------------------------
# unpack the data that has just been received;             
#---------------------------------------------------------------------
       if (timeron) timer_start(t_rdis2) end

     @info "$clusterid/$node: z=$z RECV FROM EAST --- src=$(successor[1]) -- tag = $EAST -- size=$(b_size[1])"
     @info "$clusterid/$node: z=$z RECV FROM WEST --- src=$(predecessor[1]) -- tag = $WEST -- size=$(b_size[2])"
     @info "$clusterid/$node: z=$z RECV FROM SOUTH --- src=$(successor[2]) -- tag = $NORTH -- size=$(b_size[3])"
     @info "$clusterid/$node: z=$z RECV FROM NORTH --- src=$(predecessor[2]) -- tag = $SOUTH -- size=$(b_size[4])"

       requests[1] = MPI.Irecv!(view(in_buffer, sr[1]:sr[1]+b_size[1]-1), comm_exch; source = successor[1], tag = EAST)
       requests[2] = MPI.Irecv!(view(in_buffer, sr[2]:sr[2]+b_size[2]-1), comm_exch; source = predecessor[1], tag = WEST)
       requests[3] = MPI.Irecv!(view(in_buffer, sr[3]:sr[3]+b_size[3]-1), comm_exch; source = successor[2], tag = NORTH)
       requests[4] = MPI.Irecv!(view(in_buffer, sr[4]:sr[4]+b_size[4]-1), comm_exch; source = predecessor[2], tag = SOUTH)

       MPI.Waitall(requests)

       if z==1 && clusterid==0 && node==0
        for i = eachindex(in_buffer)
       @info "$clusterid/$node: z=$z in_buffer[$i] = $(in_buffer[i])"
        end
       end

       p0 = 0
       p1 = 0
       p2 = 0
       p3 = 0

       for c = 1:ncells

            if cell_coord[1, c] == 1
               zone_cluster = zone_proc_id[iz_west[proc_zone_id[z]]]
               if zone_cluster != clusterid
                    u_face = remotecall(collect_face, 1, z, cell_low[2, c] + cell_start[2, c] + 1, cell_high[2, c] - cell_end[2,c] + 1, 
                                                            cell_low[3, c] + cell_start[3, c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                                         Val(FROM_WEST); role=:worker)

                    view(u, 0, cell_start[2,c]:cell_size[2, c]-cell_end[2,c]-1, 
                               cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, 1:5, c) .= fetch(u_face; role=:worker)
               else
                    #@info "$clusterid/$node: zone=$(proc_zone_id[z]) other_zone=$(iz_west[proc_zone_id[z]]) zone_cluster=$zone_cluster FROM_WEST"
                    for m = 1:5
                        for k = cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1
                            for j = cell_start[2,c]:cell_size[2, c]-cell_end[2,c]-1
                                i = 0
                                u[i, j, k, m, c] = in_buffer[sr[2] + p0]
                              @info "$clusterid/$node: z=$z recv-u[$i,$j,$k,$m,$c]=$(u[i,j,k,m,c]) - FROM WEST p0=$p0 sr[2]=$(sr[2])"
                                p0 = p0 + 1
                            end
                        end
                    end
               end
            end

            if cell_coord[1, c] == ncells
               zone_cluster = zone_proc_id[iz_east[proc_zone_id[z]]]
               if zone_cluster != clusterid
                    u_face = remotecall(collect_face, 1, z, cell_low[2, c] + cell_start[2,c] + 1, cell_high[2, c] - cell_end[2, c] + 1, 
                                                            cell_low[3, c] + cell_start[3,c] + 1, cell_high[3, c] - cell_end[3, c] + 1, 
                                                         Val(FROM_EAST); role=:worker)
                    view(u, cell_size[1, c]-1, cell_start[2,c]:cell_size[2, c]-cell_end[2,c]-1, 
                                               cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, 1:5, c) .= fetch(u_face; role=:worker)
               else
                    #@info "$clusterid/$node: zone=$(proc_zone_id[z]) other_zone=$(iz_east[proc_zone_id[z]]) zone_cluster=$zone_cluster FROM_EAST"
                    for m = 1:5
                        for k = cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1
                            for j = cell_start[2,c]:cell_size[2, c]-cell_end[2,c]-1
                                i = cell_size[1, c]-1
                                u[i, j, k, m, c] = in_buffer[sr[1] + p1]
                              @info "$clusterid/$node: z=$z recv-u[$i,$j,$k,$m,$c]=$(u[i,j,k,m,c]) - FROM EAST p1=$p1 sr[1]=$(sr[1])"
                                p1 = p1 + 1
                         end
                      end
                    end
               end
            end

            if cell_coord[2, c] == 1
               zone_cluster = zone_proc_id[iz_north[proc_zone_id[z]]]
               if zone_cluster != clusterid
                    u_face = remotecall(collect_face, 1, z, cell_low[1, c] + cell_start[1, c] + 1, cell_high[1, c] - cell_end[1,c] + 1, 
                                                            cell_low[3, c] + cell_start[3, c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                                         Val(FROM_NORTH); role=:worker)

                    view(u, cell_start[1,c]:cell_size[1, c]-cell_end[1,c]-1, 0, 
                            cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, 1:5, c) .= fetch(u_face; role=:worker)  
               else
                    #@info "$clusterid/$node: zone=$(proc_zone_id[z]) other_zone=$(iz_north[proc_zone_id[z]]) zone_cluster=$zone_cluster FROM_NORTH"
                    for m = 1:5
                        for k = cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1
                            j = 0
                            for i = cell_start[1,c]:cell_size[1, c]-cell_end[1,c]-1
                                u[i, j, k, m, c] = in_buffer[sr[4] + p2]
                                x = in_buffer[961 + p2]
                              @info "$clusterid/$node: z=$z recv-u[$i,$j,$k,$m,$c]=$(u[i,j,k,m,c]) - FROM NORTH p2=$p2 sr[4]=$(sr[4]) x=$x"
                                p2 = p2 + 1
                            end
                        end
                    end
               end
            end

            if cell_coord[2, c] == ncells
               zone_cluster = zone_proc_id[iz_south[proc_zone_id[z]]]
               if zone_cluster != clusterid
                  u_face = remotecall(collect_face, 1, z, cell_low[1, c] + cell_start[1,c] + 1, cell_high[1, c] - cell_end[1,c] + 1, 
                                                          cell_low[3, c] + cell_start[3,c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                                      Val(FROM_SOUTH); role=:worker)

                  view(u, cell_start[1,c]:cell_size[1,c]-cell_end[1,c]-1, cell_size[2, c]-1, 
                          cell_start[3,c]:cell_size[3,c]-cell_end[3,c]-1, 1:5, c) .= fetch(u_face; role=:worker)
               else
                  #@info "$clusterid/$node: zone=$(proc_zone_id[z]) other_zone=$(iz_south[proc_zone_id[z]]) zone_cluster=$zone_cluster FROM_SOUTH"
                  for m = 1:5
                      for k = cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1
                          j = cell_size[2, c]-1
                          for i = cell_start[1,c]:cell_size[1, c]-cell_end[1,c]-1
                              u[i, j, k, m, c] = in_buffer[sr[3] + p3]
                            @info "$clusterid/$node: z=$z recv-u[$i,$j,$k,$m,$c]=$(u[i,j,k,m,c]) - FROM SOUTH p3=$p3 sr[3]=$(sr[3])"
                              p3 = p3 + 1
                          end
                      end
                  end
               end
            end

#---------------------------------------------------------------------
#      cells loop
#---------------------------------------------------------------------
       end
       if (timeron) timer_stop(t_rdis2) end


        return nothing
end

   