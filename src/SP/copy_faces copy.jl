#---------------------------------------------------------------------
#---------------------------------------------------------------------


function copy_faces()
         
#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# this function copies the face values of a variable defined on a set 
# of cells to the overlap locations of the adjacent sets of cells. 
# Because a set of cells interfaces in each direction with exactly one 
# other set, we only need to fill six different buffers. We could try to 
# overlap communication with computation, by computing
# some internal values while communicating boundary values, but this
# adds so much overhead that it's not clearly useful. 
#---------------------------------------------------------------------

#       use sp_data
#       use mpinpb

#       implicit none

#       integer i, j, k, c, m, requests[0:11], p0, p1,  
#               p2, p3, p4, p5, b_size[0:5], ss[0:5],  
#               sr[0:5], ERROR, statuses[MPI_STATUS_SIZE, 0:11]

       requests = Array{MPI.Request}(undef,12)

#---------------------------------------------------------------------
#      exit immediately if there are no faces to be copied           
#---------------------------------------------------------------------
       if no_nodes == 1
          compute_rhs()
          return nothing
       end

#---------------------------------------------------------------------
# because the difference stencil for the diagonalized scheme is 
# orthogonal, we do not have to perform the staged copying of faces, 
# but can send all face information simultaneously to the neighboring 
# cells in all directions          
#---------------------------------------------------------------------
       if (timeron) timer_start(t_bpack) end
       p0 = 1
       p1 = 1
       p2 = 1
       p3 = 1
       p4 = 1
       p5 = 1

       for c = 1:ncells
          for m = 1:5

#---------------------------------------------------------------------
#            fill the buffer to be sent to eastern neighbors (i-dir)
#---------------------------------------------------------------------
             if cell_coord[1, c] != ncells
                for k = 0:cell_size[3, c]-1
                   for j = 0:cell_size[2, c]-1
                      for i = cell_size[1, c]-2:cell_size[1, c]-1
                        out_buffer[1][p0] = u[i, j, k, m, c]
                        p0 = p0 + 1
                      end
                   end
                end
             end

#---------------------------------------------------------------------
#            fill the buffer to be sent to western neighbors 
#---------------------------------------------------------------------
             if cell_coord[1, c] != 1
                for k = 0:cell_size[3, c]-1
                   for j = 0:cell_size[2, c]-1
                      for i = 0:1
                         out_buffer[2][p1] = u[i, j, k, m, c]
                         p1 = p1 + 1
                      end
                   end
                end


             end

#---------------------------------------------------------------------
#            fill the buffer to be sent to northern neighbors (j_dir)
#---------------------------------------------------------------------
             if cell_coord[2, c] != ncells
                for k = 0:cell_size[3, c]-1
                   for j = cell_size[2, c]-2:cell_size[2, c]-1
                      for i = 0:cell_size[1, c]-1
                         out_buffer[3][p2] = u[i, j, k, m, c]
                         p2 = p2 + 1
                      end
                   end
                end
             end

#---------------------------------------------------------------------
#            fill the buffer to be sent to southern neighbors 
#---------------------------------------------------------------------
             if cell_coord[2, c] != 1
                for k = 0:cell_size[3, c]-1
                   for j = 0:1
                      for i = 0:cell_size[1, c]-1
                         out_buffer[4][p3] = u[i, j, k, m, c]
                         p3 = p3 + 1
                      end
                   end
                end
             end

#---------------------------------------------------------------------
#            fill the buffer to be sent to top neighbors (k-dir)
#---------------------------------------------------------------------
             if cell_coord[3, c] != ncells
                for k = cell_size[3, c]-2:cell_size[3, c]-1
                   for j = 0:cell_size[2, c]-1
                      for i = 0:cell_size[1, c]-1
                         out_buffer[5][p4] = u[i, j, k, m, c]
                         p4 = p4 + 1
                      end
                   end
                end
             end

#---------------------------------------------------------------------
#            fill the buffer to be sent to bottom neighbors
#---------------------------------------------------------------------
             if cell_coord[3, c] != 1
                 for k = 0:1
                    for j = 0:cell_size[2, c]-1
                       for i = 0:cell_size[1, c]-1
                          out_buffer[6][p5] = u[i, j, k, m, c]
                          p5 = p5 + 1
                       end
                    end
                 end
              end

#---------------------------------------------------------------------
#          m loop
#---------------------------------------------------------------------
           end

#---------------------------------------------------------------------
#       cell loop
#---------------------------------------------------------------------
        end
       if (timeron) timer_stop(t_bpack) end

       if (timeron) timer_start(t_exch) end

       requests[1] = MPI.Irecv!(in_buffer[1], comm_rhs; source = sucessor[1], tag = WEST)
       requests[2] = MPI.Irecv!(in_buffer[2], comm_rhs; source = predecessor[1], tag = EAST)
       requests[3] = MPI.Irecv!(in_buffer[3], comm_rhs; source = sucessor[2], tag = SOUTH)
       requests[4] = MPI.Irecv!(in_buffer[4], comm_rhs; source = predecessor[2], tag = NORTH)
       requests[5] = MPI.Irecv!(in_buffer[5], comm_rhs; source = sucessor[3], tag = BOTTOM)
       requests[6] = MPI.Irecv!(in_buffer[6], comm_rhs; source = predecessor[3], tag = TOP)

       requests[7] = MPI.Isend(out_buffer[1], comm_rhs; dest = sucessor[1], tag = EAST)
       requests[8] = MPI.Isend(out_buffer[2], comm_rhs; dest = predecessor[1], tag = WEST)
       requests[9] = MPI.Isend(out_buffer[3], comm_rhs; dest = sucessor[2], tag = NORTH)
       requests[10] = MPI.Isend(out_buffer[4], comm_rhs; dest = predecessor[2], tag = SOUTH)
       requests[11] = MPI.Isend(out_buffer[5], comm_rhs; dest = sucessor[3], tag = TOP)
       requests[12] = MPI.Isend(out_buffer[6], comm_rhs; dest = predecessor[3], tag = BOTTOM)

       MPI.Waitall(requests)

       if (timeron) timer_stop(t_exch) end

#---------------------------------------------------------------------
# unpack the data that has just been received;             
#---------------------------------------------------------------------
       if (timeron) timer_start(t_bpack) end
       p0 = 1
       p1 = 1
       p2 = 1
       p3 = 1
       p4 = 1
       p5 = 1

       for c = 1:ncells
          for m = 1:5

             if cell_coord[1, c] != 1
                for k = 0:cell_size[3, c]-1
                   for j = 0:cell_size[2, c]-1
                      for i = -2:-1
                         u[i, j, k, m, c] = in_buffer[1][p0]
                         p0 = p0 + 1
                      end
                   end
                end
             end

             if cell_coord[1, c] != ncells
                for k = 0:cell_size[3, c]-1
                   for j = 0:cell_size[2, c]-1
                      for i = cell_size[1, c]:cell_size[1, c]+1
                         u[i, j, k, m, c] = in_buffer[2][p1]
                         p1 = p1 + 1
                      end
                   end
                end
             end

             if cell_coord[2, c] != 1
                for k = 0:cell_size[3, c]-1
                   for j = -2:-1
                      for i = 0:cell_size[1, c]-1
                         u[i, j, k, m, c] = in_buffer[3][p2]
                         p2 = p2 + 1
                      end
                   end
                end

             end

             if cell_coord[2, c] != ncells
                for k = 0:cell_size[3, c]-1
                   for j = cell_size[2, c]:cell_size[2, c]+1
                      for i = 0:cell_size[1, c]-1
                         u[i, j, k, m, c] = in_buffer[4][p3]
                         p3 = p3 + 1
                      end
                   end
                end
             end

             if cell_coord[3, c] != 1
                for k = -2:-1
                   for j = 0:cell_size[2, c]-1
                      for i = 0:cell_size[1, c]-1
                         u[i, j, k, m, c] = in_buffer[5][p4]
                         p4 = p4 + 1
                      end
                   end
                end
             end

             if cell_coord[3, c] != ncells
                for k = cell_size[3, c]:cell_size[3, c]+1
                   for j = 0:cell_size[2, c]-1
                      for i = 0:cell_size[1, c]-1
                         u[i, j, k, m, c] = in_buffer[6][p5]
                         p5 = p5 + 1
                      end
                   end
                end
             end

#---------------------------------------------------------------------
#         m loop            
#---------------------------------------------------------------------
          end

#---------------------------------------------------------------------
#      cells loop
#---------------------------------------------------------------------
       end
       if (timeron) timer_stop(t_bpack) end

#---------------------------------------------------------------------
# now that we have all the data, compute the rhs
#---------------------------------------------------------------------
       compute_rhs()

return nothing
end
