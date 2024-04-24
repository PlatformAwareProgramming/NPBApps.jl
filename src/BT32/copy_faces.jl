
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

function copy_faces(z, sr, 
                     b_size,
                     cell_coord,
                     cell_size,
                     cell_start,
                     cell_end,
                     forcing,        
                     u,
                     rhs,
                     in_buffer,
                     out_buffer,
                     us,
                     vs,
                     ws,
                     qs,
                     rho_i,
                     square,
                     timeron,
                     dt,
                     ::Val{ncells},
                     tx2,
                     ty2,
                     tz2,
                     dx1tx1,
                     dx2tx1,
                     dx3tx1,
                     dx4tx1,
                     dx5tx1,
                     dy1ty1,
                     dy2ty1,
                     dy3ty1,
                     dy4ty1,
                     dy5ty1,
                     dz1tz1,
                     dz2tz1,
                     dz3tz1,
                     dz4tz1,
                     dz5tz1,
                     xxcon2,
                     xxcon3,
                     xxcon4,
                     xxcon5,
                     yycon2,
                     yycon3,
                     yycon4,
                     yycon5,
                     zzcon2,
                     zzcon3,
                     zzcon4,
                     zzcon5,
                     _::Val{1}, 
                     comm_rhs,
                     predecessor,
                     successor,
                     requests,
                     ) where ncells

         compute_rhs(cell_size,
                     cell_start,
                     cell_end,
                     forcing,           
                     u,
                     rhs,
                     us,
                     vs,
                     ws,
                     qs,
                     rho_i,
                     square,
                     timeron,
                     dt,
                     ncells,
                     tx2,
                     ty2,
                     tz2,
                     dx1tx1,
                     dx2tx1,
                     dx3tx1,
                     dx4tx1,
                     dx5tx1,
                     dy1ty1,
                     dy2ty1,
                     dy3ty1,
                     dy4ty1,
                     dy5ty1,
                     dz1tz1,
                     dz2tz1,
                     dz3tz1,
                     dz4tz1,
                     dz5tz1,
                     xxcon2,
                     xxcon3,
                     xxcon4,
                     xxcon5,
                     yycon2,
                     yycon3,
                     yycon4,
                     yycon5,
                     zzcon2,
                     zzcon3,
                     zzcon4,
                     zzcon5,
      )

end



function copy_faces(z, ss, 
                     sr, 
                     b_size,
                     cell_coord,
                     cell_size,
                     cell_start,
                     cell_end,
                     forcing,        
                     u,
                     rhs,
                     in_buffer,
                     out_buffer,
                     us,
                     vs,
                     ws,
                     qs,
                     rho_i,
                     square,
                     timeron,
                     dt,
                     ::Val{ncells},
                     tx2,
                     ty2,
                     tz2,
                     dx1tx1,
                     dx2tx1,
                     dx3tx1,
                     dx4tx1,
                     dx5tx1,
                     dy1ty1,
                     dy2ty1,
                     dy3ty1,
                     dy4ty1,
                     dy5ty1,
                     dz1tz1,
                     dz2tz1,
                     dz3tz1,
                     dz4tz1,
                     dz5tz1,
                     xxcon2,
                     xxcon3,
                     xxcon4,
                     xxcon5,
                     yycon2,
                     yycon3,
                     yycon4,
                     yycon5,
                     zzcon2,
                     zzcon3,
                     zzcon4,
                     zzcon5,
                     _::Val{no_nodes}, 
                     comm_rhs,
                     predecessor,
                     successor,
                     requests,
                     ) where {no_nodes, ncells}


##@info "$clusterid/$node: copy faces 0 - z=$z" 
                     
#---------------------------------------------------------------------
#     because the difference stencil for the diagonalized scheme is 
#     orthogonal, we do not have to perform the staged copying of faces, 
#     but can send all face information simultaneously to the neighboring 
#     cells in all directions          
#---------------------------------------------------------------------
      if (timeron) timer_start(t_bpack) end
      p0 = 0
      p1 = 0
      p2 = 0
      p3 = 0
      p4 = 0
      p5 = 0

      for c = 1:ncells
         ##@info "$clusterid/$node: copy faces 1 BEGIN - c=$c" 

#---------------------------------------------------------------------
#     fill the buffer to be sent to eastern neighbors (i-dir)
#---------------------------------------------------------------------
         if cell_coord[1, c] != ncells
            ##@info "$clusterid/$node: copy faces inner 2.1 - z=$z c=$c" 
            for k = 0:cell_size[3, c]-1
               for j = 0:cell_size[2, c]-1
                  for i = cell_size[1, c]-2:cell_size[1, c]-1
                     for m = 1:5
                        out_buffer[ss[1]+p0] = u[m, i, j, k, c]
                        p0 = p0 + 1
                     end
                  end
               end
            end
            ##@info "$clusterid/$node: copy faces inner 2.2 - z=$z c=$c" 
         end

#---------------------------------------------------------------------
#     fill the buffer to be sent to western neighbors 
#---------------------------------------------------------------------
         if cell_coord[1, c] != 1
            ##@info "$clusterid/$node: copy faces inner 3.1 - z=$z c=$c" 
            for k = 0:cell_size[3, c]-1
               for j = 0:cell_size[2, c]-1
                  for i = 0:1
                     for m = 1:5
                        out_buffer[ss[2]+p1] = u[m, i, j, k, c]
                        p1 = p1 + 1
                     end
                  end
               end
            end
            ##@info "$clusterid/$node: copy faces inner 3.2 - z=$z c=$c" 
         end

#---------------------------------------------------------------------
#     fill the buffer to be sent to northern neighbors (j_dir)
#---------------------------------------------------------------------
         if cell_coord[2, c] != ncells
            ##@info "$clusterid/$node: copy faces inner 4.1 - z=$z c=$c" 
            for k = 0:cell_size[3, c]-1
               for j = cell_size[2, c]-2:cell_size[2, c]-1
                  for i = 0:cell_size[1, c]-1
                     for m = 1:5
                        out_buffer[ss[3]+p2] = u[m, i, j, k, c]
                        p2 = p2 + 1
                     end
                  end
               end
            end
            ##@info "$clusterid/$node: copy faces inner 4.2 - z=$z c=$c" 
         end

#---------------------------------------------------------------------
#     fill the buffer to be sent to southern neighbors 
#---------------------------------------------------------------------
         if cell_coord[2, c] != 1
            ##@info "$clusterid/$node: copy faces inner 5.1 - z=$z c=$c" 
            for k = 0:cell_size[3, c]-1
               for j = 0:1
                  for i = 0:cell_size[1, c]-1
                     for m = 1:5
                        out_buffer[ss[4]+p3] = u[m, i, j, k, c]
                        p3 = p3 + 1
                     end
                  end
               end
            end
            ##@info "$clusterid/$node: copy faces inner 5.2 - z=$z c=$c" 
         end

#---------------------------------------------------------------------
#     fill the buffer to be sent to top neighbors (k-dir)
#---------------------------------------------------------------------
         if cell_coord[3, c] != ncells
            ##@info "$clusterid/$node: copy faces inner 101.1 - z=$z c=$c" 
            for k = cell_size[3, c]-2:cell_size[3, c]-1
               for j = 0:cell_size[2, c]-1
                  for i = 0:cell_size[1, c]-1
                     for m = 1:5
                        # ##@info "$clusterid/$node: zone=$z length(out_buffer)=$(length(out_buffer)) out_buffer[$(ss[5]+p4)] --- m=$m i=$i j=$j k=$k ss[5]=$(ss[5]) p4=$p4 - BEGIN"
                        out_buffer[ss[5]+p4] = u[m, i, j, k, c]
                        # ##@info "$clusterid/$node: zone=$z length(out_buffer)=$(length(out_buffer)) out_buffer[$(ss[5]+p4)] --- m=$m i=$i j=$j k=$k ss[5]=$(ss[5]) p4=$p4 - END"
                        p4 = p4 + 1
                     end
                  end
               end
            end
            ##@info "$clusterid/$node: copy faces inner 101.2 - z=$z c=$c" 
         end

#---------------------------------------------------------------------
#     fill the buffer to be sent to bottom neighbors
#---------------------------------------------------------------------
         if cell_coord[3, c] != 1
            ##@info "$clusterid/$node: copy faces inner 102.1 - z=$z c=$c" 
            for k = 0:1
               for j = 0:cell_size[2, c]-1
                  for i = 0:cell_size[1, c]-1
                     for m = 1:5
                        out_buffer[ss[6]+p5] = u[m, i, j, k, c]
                        p5 = p5 + 1
                     end
                  end
               end
            end
            ##@info "$clusterid/$node: copy faces inner 102.2 - z=$z c=$c" 
         end
         ##@info "$clusterid/$node: copy faces 1 END - c=$c" 
#---------------------------------------------------------------------
#     cell loop
#---------------------------------------------------------------------
      end
      if (timeron) timer_stop(t_bpack) end

      if (timeron) timer_start(t_exch) end

      ##@info "$clusterid/$node: copy faces 6 - z=$z" 

      requests[1] = MPI.Irecv!(view(in_buffer, sr[1]:sr[1]+b_size[1]-1), comm_rhs; source = successor[1], tag = WEST)
      requests[2] = MPI.Irecv!(view(in_buffer, sr[2]:sr[2]+b_size[2]-1), comm_rhs; source = predecessor[1], tag = EAST)
      requests[3] = MPI.Irecv!(view(in_buffer, sr[3]:sr[3]+b_size[3]-1), comm_rhs; source = successor[2], tag = SOUTH)
      requests[4] = MPI.Irecv!(view(in_buffer, sr[4]:sr[4]+b_size[4]-1), comm_rhs; source = predecessor[2], tag = NORTH)
      requests[5] = MPI.Irecv!(view(in_buffer, sr[5]:sr[5]+b_size[5]-1), comm_rhs; source = successor[3], tag = BOTTOM)
      requests[6] = MPI.Irecv!(view(in_buffer, sr[6]:sr[6]+b_size[6]-1), comm_rhs; source = predecessor[3], tag = TOP)

      requests[7] = MPI.Isend(view(out_buffer, ss[1]:ss[1]+b_size[1]-1), comm_rhs; dest = successor[1], tag = EAST)
      requests[8] = MPI.Isend(view(out_buffer, ss[2]:ss[2]+b_size[2]-1), comm_rhs; dest = predecessor[1], tag = WEST)
      requests[9] = MPI.Isend(view(out_buffer, ss[3]:ss[3]+b_size[3]-1), comm_rhs; dest = successor[2], tag = NORTH)
      requests[10] = MPI.Isend(view(out_buffer, ss[4]:ss[4]+b_size[4]-1), comm_rhs; dest = predecessor[2], tag = SOUTH)
      requests[11] = MPI.Isend(view(out_buffer, ss[5]:ss[5]+b_size[5]-1), comm_rhs; dest = successor[3], tag = TOP)
      requests[12] = MPI.Isend(view(out_buffer, ss[6]:ss[6]+b_size[6]-1), comm_rhs; dest = predecessor[3], tag = BOTTOM)

      ##@info "$clusterid/$node: copy faces 7 - z=$z" 

      #try
         #@info "$clusterid/$node: z=$z ARRIVE WAIT ALL"
      #   lock(lock_comm)
      #   #@info "$clusterid/$node: z=$z ENTER WAIT ALL"
         MPI.Waitall(requests)
         #@info "$clusterid/$node: z=$z DEPART WAIT ALL"
      #finally
      #   unlock(lock_comm)
      #end

      ##@info "$clusterid/$node: copy faces 8 - z=$z" 

      if (timeron) timer_stop(t_exch) end

#---------------------------------------------------------------------
#     unpack the data that has just been received;             
#---------------------------------------------------------------------
      if (timeron) timer_start(t_bpack) end
      p0 = 0
      p1 = 0
      p2 = 0
      p3 = 0
      p4 = 0
      p5 = 0

      for c = 1:ncells

         if cell_coord[1, c] != 1
            for k = 0:cell_size[3, c]-1
               for j = 0:cell_size[2, c]-1
                  for i = -2:-1
                     for m = 1:5
                        u[m, i, j, k, c] = in_buffer[sr[2]+p0]
                        p0 = p0 + 1
                     end
                  end
               end
            end
         end

         if cell_coord[1, c] != ncells
            for k = 0:cell_size[3, c]-1
               for j = 0:cell_size[2, c]-1
                  for i = cell_size[1, c]:cell_size[1, c]+1
                     for m = 1:5
                        u[m, i, j, k, c] = in_buffer[sr[1]+p1]
                        p1 = p1 + 1
                     end
                  end
               end
            end
         end

         if cell_coord[2, c] != 1
            for k = 0:cell_size[3, c]-1
               for j = -2:-1
                  for i = 0:cell_size[1, c]-1
                     for m = 1:5
                        u[m, i, j, k, c] = in_buffer[sr[4]+p2]
                        p2 = p2 + 1
                     end
                  end
               end
            end
         end

         if cell_coord[2, c] != ncells
            for k = 0:cell_size[3, c]-1
               for j = cell_size[2, c]:cell_size[2, c]+1
                  for i = 0:cell_size[1, c]-1
                     for m = 1:5
                        u[m, i, j, k, c] = in_buffer[sr[3]+p3]
                        p3 = p3 + 1
                     end
                  end
               end
            end
         end

         if cell_coord[3, c] != 1
            for k = -2:-1
               for j = 0:cell_size[2, c]-1
                  for i = 0:cell_size[1, c]-1
                     for m = 1:5
                        u[m, i, j, k, c] = in_buffer[sr[6]+p4]
                        p4 = p4 + 1
                     end
                  end
               end
            end
         end

         if cell_coord[3, c] != ncells
            for k = cell_size[3, c]:cell_size[3, c]+1
               for j = 0:cell_size[2, c]-1
                  for i = 0:cell_size[1, c]-1
                     for m = 1:5
                        u[m, i, j, k, c] = in_buffer[sr[5]+p5]
                        p5 = p5 + 1
                     end
                  end
               end
            end
         end

#---------------------------------------------------------------------
#     cells loop
#---------------------------------------------------------------------
      end
      if (timeron) timer_stop(t_bpack) end

      ##@info "$clusterid/$node: copy faces 13 - z=$z" 

#---------------------------------------------------------------------
#     do the rest of the rhs that uses the copied face values          
#---------------------------------------------------------------------
      compute_rhs(cell_size,
                  cell_start,
                  cell_end,
                  forcing,           
                  u,
                  rhs,
                  us,
                  vs,
                  ws,
                  qs,
                  rho_i,
                  square,
                  timeron,
                  dt,
                  ncells,
                  tx2,
                  ty2,
                  tz2,
                  dx1tx1,
                  dx2tx1,
                  dx3tx1,
                  dx4tx1,
                  dx5tx1,
                  dy1ty1,
                  dy2ty1,
                  dy3ty1,
                  dy4ty1,
                  dy5ty1,
                  dz1tz1,
                  dz2tz1,
                  dz3tz1,
                  dz4tz1,
                  dz5tz1,
                  xxcon2,
                  xxcon3,
                  xxcon4,
                  xxcon5,
                  yycon2,
                  yycon3,
                  yycon4,
                  yycon5,
                  zzcon2,
                  zzcon3,
                  zzcon4,
                  zzcon5,
               )

      ##@info "$clusterid/$node: copy faces 14 - z=$z" 

      return nothing
end
