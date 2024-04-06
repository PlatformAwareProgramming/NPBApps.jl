function exch_qbc(_::Val{ncells},
                    proc_num_zones,
                    cell_coord,
                    cell_size,
                    cell_start,
                    cell_end,
                    cell_low,
                    cell_high,
                    u,
                    timeron,
                    ) where {ncells}  

    for iz = 1:proc_num_zones
        exch_qbc_send(Val(ncells),
                        iz,
                        cell_coord[iz],
                        cell_size[iz],
                        cell_start[iz],
                        cell_end[iz],
                        cell_low[iz],
                        cell_high[iz],
                        u[iz],
                        timeron,)
    end
 
    for iz = 1:proc_num_zones
        exch_qbc_recv(Val(ncells),
                        iz,
                        cell_coord[iz],
                        cell_size[iz],
                        cell_start[iz],
                        cell_end[iz],
                        cell_low[iz],
                        cell_high[iz],
                        u[iz],
                        timeron,)
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


function exch_qbc_send(_::Val{ncells},
                            z,
                            cell_coord,
                            cell_size,
                            cell_start,
                            cell_end,
                            cell_low,
                            cell_high,
                            u,
                            timeron,
                           ) where {ncells}  
                           
#---------------------------------------------------------------------
# because the difference stencil for the diagonalized scheme is 
# orthogonal, we do not have to perform the staged copying of faces, 
# but can send all face information simultaneously to the neighboring 
# cells in all directions          
#---------------------------------------------------------------------
       if (timeron) timer_start(t_rdis1) end
 
       for c = 1:ncells

#---------------------------------------------------------------------
#            fill the buffer to be sent to eastern neighbors (i-dir)
#---------------------------------------------------------------------
            if cell_coord[1, c] == ncells
               # outer face (inter-zone)
               u_face = view(u, cell_size[1, c]-2, cell_start[2,c]:cell_size[2, c]-cell_end[2,c]-1, 
                                cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, 1:5, c)
               remotecall(deposit_face, 1, z, cell_low[2, c] + cell_start[2,c] + 1, cell_high[2, c] - cell_end[2,c] + 1, 
                                                 cell_low[3, c] + cell_start[3,c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                           u_face, Val(1); role=:worker)
            end

#---------------------------------------------------------------------
#            fill the buffer to be sent to western neighbors 
#---------------------------------------------------------------------
             if cell_coord[1, c] == 1
                # outer face (inter-zone)
               u_face = view(u, 1, cell_start[2,c]:cell_size[2, c]-cell_end[2,c]-1, 
                                   cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, 1:5, c)
               remotecall(deposit_face, 1, z, cell_low[2, c] + cell_start[2,c] + 1, cell_high[2, c] - cell_end[2,c] + 1, 
                                                 cell_low[3, c] + cell_start[3,c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                           u_face, Val(2); role=:worker)
             end

#---------------------------------------------------------------------
#            fill the buffer to be sent to northern neighbors (j_dir)
#---------------------------------------------------------------------
             if cell_coord[2, c] == ncells
                # outer face (inter-zone)
                u_face = view(u, cell_start[1,c]:cell_size[1, c]-cell_end[1,c]-1, cell_size[2, c]-2, 
                                 cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, 1:5, c)
                remotecall(deposit_face, 1, z, cell_low[1, c] + cell_start[1,c] + 1, cell_high[1, c] - cell_end[1,c] + 1, 
                                                  cell_low[3, c] + cell_start[3,c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                            u_face, Val(3); role=:worker)
             end

#---------------------------------------------------------------------
#            fill the buffer to be sent to southern neighbors 
#---------------------------------------------------------------------
             if cell_coord[2, c] == 1
                u_face = view(u, cell_start[1,c]:cell_size[1, c]-cell_end[1,c]-1, 1, 
                                 cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, 1:5, c)
                remotecall(deposit_face, 1, z, cell_low[1, c] + cell_start[1,c] + 1, cell_high[1, c] - cell_end[1,c] + 1, 
                                                  cell_low[3, c] + cell_start[3,c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                            u_face, Val(4); role=:worker)
            end

#---------------------------------------------------------------------
#       cell loop
#---------------------------------------------------------------------
       end

       if (timeron) timer_stop(t_rdis1) end

end


function exch_qbc_recv(_::Val{ncells},
                            z,
                            cell_coord,
                            cell_size,
                            cell_start,
                            cell_end,
                            cell_low,
                            cell_high,
                            u,
                            timeron,
                     ) where {ncells}  

#---------------------------------------------------------------------
# unpack the data that has just been received;             
#---------------------------------------------------------------------
       if (timeron) timer_start(t_rdis2) end

       for c = 1:ncells

             if cell_coord[1, c] == 1
                u_face = remotecall(collect_face, 1, z, cell_low[2, c] + cell_start[2,c] + 1, cell_high[2, c] - cell_end[2, c] + 1, 
                                                           cell_low[3, c] + cell_start[3,c] + 1, cell_high[3, c] - cell_end[3, c] + 1, 
                                                     Val(2); role=:worker)
                view(u, 0, cell_start[2,c]:cell_size[2, c]-cell_end[2,c]-1, 
                              cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, 1:5, c) .= fetch(u_face; role=:worker)
             end

             if cell_coord[1, c] == ncells
                u_face = remotecall(collect_face, 1, z, cell_low[2, c] + cell_start[2, c] + 1, cell_high[2, c] - cell_end[2,c] + 1, 
                                                           cell_low[3, c] + cell_start[3, c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                                     Val(1); role=:worker)
                view(u, cell_size[1, c]-1, cell_start[2,c]:cell_size[2, c]-cell_end[2,c]-1, 
                                                 cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, 1:5, c) .= fetch(u_face; role=:worker)
             end

             if cell_coord[2, c] == 1
               u_face = remotecall(collect_face, 1, z, cell_low[1, c] + cell_start[1, c] + 1, cell_high[1, c] - cell_end[1,c] + 1, 
                                                          cell_low[3, c] + cell_start[3, c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                                    Val(4); role=:worker)
               view(u, cell_start[1,c]:cell_size[1, c]-cell_end[1,c]-1, 0, 
                          cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, 1:5, c) .= fetch(u_face; role=:worker)  
            end

             if cell_coord[2, c] == ncells
                u_face = remotecall(collect_face, 1, z, cell_low[1, c] + cell_start[1,c] + 1, cell_high[1, c] - cell_end[1,c] + 1, 
                                                           cell_low[3, c] + cell_start[3,c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                                     Val(3); role=:worker)
                view(u, cell_start[1,c]:cell_size[1,c]-cell_end[1,c]-1, cell_size[2, c]-1, 
                           cell_start[3,c]:cell_size[3,c]-cell_end[3,c]-1, 1:5, c) .= fetch(u_face; role=:worker)
              end

#---------------------------------------------------------------------
#      cells loop
#---------------------------------------------------------------------
       end
       if (timeron) timer_stop(t_rdis2) end


        return nothing
end
