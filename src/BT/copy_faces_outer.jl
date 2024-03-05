
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


function copy_faces_outer_1( z, 
                     cell_coord,
                     cell_size,
                     cell_start,
                     cell_end,
                     cell_low,
                     cell_high,
                     u,
                     timeron,
                     ::Val{ncells}
                     ) where {ncells}


      # @info "$clusterid/$node: copy faces 0 - z=$z" 
                     
#---------------------------------------------------------------------
#     because the difference stencil for the diagonalized scheme is 
#     orthogonal, we do not have to perform the staged copying of faces, 
#     but can send all face information simultaneously to the neighboring 
#     cells in all directions          
#---------------------------------------------------------------------
      if (timeron) timer_start(t_rdis1) end

      for c = 1:ncells
         # @info "$clusterid/$node: copy faces 1 BEGIN - c=$c" 

         @info "$clusterid/$node: copy faces 1.0 - z=$z c=$c ncells=$ncells --- cell_coord[1, $c] == $(cell_coord[1,c]) / cell_coord[2, $c] == $(cell_coord[2,c])"  

#---------------------------------------------------------------------
#     fill the buffer to be sent to eastern neighbors (i-dir)
#---------------------------------------------------------------------
         if cell_coord[1, c] == ncells
               # @info "$clusterid/$node: copy faces outer 2.1 - z=$z c=$c" 
                # outer face (inter-zone)
               u_face = view(u, 1:5, cell_size[1, c]-2, cell_start[2,c]:cell_size[2, c]-cell_end[2,c]-1, 
                                                        cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, c)
               remotecall(deposit_face, 1, z, cell_low[2, c] + cell_start[2,c] + 1, cell_high[2, c] - cell_end[2,c] + 1, 
                                                cell_low[3, c] + cell_start[3,c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                          u_face, Val(1); role=:worker)
               # @info "$clusterid/$node: copy faces outer 2.2 - z=$z c=$c" 
         end

#---------------------------------------------------------------------
#     fill the buffer to be sent to western neighbors 
#---------------------------------------------------------------------
         if cell_coord[1, c] == 1
            # @info "$clusterid/$node: copy faces outer 3.1 - z=$z c=$c" 
            # outer face (inter-zone)
            u_face = view(u, 1:5, 1, cell_start[2,c]:cell_size[2, c]-cell_end[2,c]-1, 
                                     cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, c)
            remotecall(deposit_face, 1, z, cell_low[2, c] + cell_start[2,c] + 1, cell_high[2, c] - cell_end[2,c] + 1, 
                                             cell_low[3, c] + cell_start[3,c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                       u_face, Val(2); role=:worker)
            # @info "$clusterid/$node: copy faces outer 3.2 - z=$z c=$c" 
      end

#---------------------------------------------------------------------
#     fill the buffer to be sent to northern neighbors (j_dir)
#---------------------------------------------------------------------
         if cell_coord[2, c] == ncells
            # @info "$clusterid/$node: copy faces outer 4.1 - z=$z c=$c" 
            # outer face (inter-zone)
            u_face = view(u, 1:5, cell_start[1,c]:cell_size[1, c]-cell_end[1,c]-1, cell_size[2, c]-2, 
                                  cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, c)
            remotecall(deposit_face, 1, z, cell_low[1, c] + cell_start[1,c] + 1, cell_high[1, c] - cell_end[1,c] + 1, 
                                             cell_low[3, c] + cell_start[3,c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                       u_face, Val(3); role=:worker)
            # @info "$clusterid/$node: copy faces outer 4.2 - z=$z c=$c" 
         end

#---------------------------------------------------------------------
#     fill the buffer to be sent to southern neighbors 
#---------------------------------------------------------------------
         if cell_coord[2, c] == 1
            # @info "$clusterid/$node: copy faces outer 5.1 - z=$z c=$c" 
            u_face = view(u, 1:5, cell_start[1,c]:cell_size[1, c]-cell_end[1,c]-1, 1, 
                                  cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, c)
            remotecall(deposit_face, 1, z, cell_low[1, c] + cell_start[1,c] + 1, cell_high[1, c] - cell_end[1,c] + 1, 
                                             cell_low[3, c] + cell_start[3,c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                       u_face, Val(4); role=:worker)                                       
            # @info "$clusterid/$node: copy faces outer 5.2 - z=$z c=$c" 
         end

         # @info "$clusterid/$node: copy faces 1 END - c=$c" 
#---------------------------------------------------------------------
#     cell loop
#---------------------------------------------------------------------
      end

      if (timeron) timer_stop(t_rdis1) end

   end

   function copy_faces_outer_2(z, 
                       cell_coord,
                       cell_size,
                       cell_start,
                       cell_end,
                       cell_low,
                       cell_high,
                       u,
                       timeron,
                       ::Val{ncells},
                     ) where {ncells}


      if (timeron) timer_start(t_rdis2) end

#---------------------------------------------------------------------
#     unpack the data that has just been received;             
#---------------------------------------------------------------------
   
      for c = 1:ncells

         if cell_coord[1, c] == 1
            # @info "$clusterid/$node: copy faces 9.1 - z=$z c=$c" 
            u_face = remotecall(collect_face, 1, z, cell_low[2, c] + cell_start[2,c] + 1, cell_high[2, c] - cell_end[2, c] + 1, 
                                                                        cell_low[3, c] + cell_start[3,c] + 1, cell_high[3, c] - cell_end[3, c] + 1, 
                                                               Val(2); role=:worker)
            view(u, 1:5, 0, cell_start[2,c]:cell_size[2, c]-cell_end[2,c]-1, 
                            cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, c) .= fetch(u_face; role=:worker)
            # @info "$clusterid/$node: copy faces 9.2 - z=$z c=$c" 
         end

         if cell_coord[1, c] == ncells
            # @info "$clusterid/$node: copy faces 10.1 - z=$z c=$c" 
            u_face = remotecall(collect_face, 1, z, cell_low[2, c] + cell_start[2, c] + 1, cell_high[2, c] - cell_end[2,c] + 1, 
                                                                     cell_low[3, c] + cell_start[3, c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                                               Val(1); role=:worker)
            view(u, 1:5, cell_size[1, c]-1, cell_start[2,c]:cell_size[2, c]-cell_end[2,c]-1, 
                                            cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, c) .= fetch(u_face; role=:worker)
            # @info "$clusterid/$node: copy faces 10.2 - z=$z c=$c" 
         end

         if cell_coord[2, c] == 1
            # @info "$clusterid/$node: copy faces 11.1 - z=$z c=$c" 
            u_face = remotecall(collect_face, 1, z, cell_low[1, c] + cell_start[1, c] + 1, cell_high[1, c] - cell_end[1,c] + 1, 
                                                                     cell_low[3, c] + cell_start[3, c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                                               Val(4); role=:worker)
            view(u, 1:5, cell_start[1,c]:cell_size[1, c]-cell_end[1,c]-1, 0, 
                         cell_start[3,c]:cell_size[3, c]-cell_end[3,c]-1, c) .= fetch(u_face; role=:worker)
            # @info "$clusterid/$node: copy faces 11.2 - z=$z c=$c" 
         end

         if cell_coord[2, c] == ncells
            # @info "$clusterid/$node: copy faces 12.1 - z=$z c=$c" 
            u_face = remotecall(collect_face, 1, z, cell_low[1, c] + cell_start[1,c] + 1, cell_high[1, c] - cell_end[1,c] + 1, 
                                                                     cell_low[3, c] + cell_start[3,c] + 1, cell_high[3, c] - cell_end[3,c] + 1, 
                                                               Val(3); role=:worker)
            view(u, 1:5, cell_start[1,c]:cell_size[1,c]-cell_end[1,c]-1, cell_size[2, c]-1, 
                         cell_start[3,c]:cell_size[3,c]-cell_end[3,c]-1, c) .= fetch(u_face; role=:worker)
            # @info "$clusterid/$node: copy faces 12.2 - z=$z c=$c" 
         end

#---------------------------------------------------------------------
#     cells loop
#---------------------------------------------------------------------
      end

      # @info "$clusterid/$node: copy faces 13 - z=$z" 

      if (timeron) timer_stop(t_rdis2) end

      return nothing
end
