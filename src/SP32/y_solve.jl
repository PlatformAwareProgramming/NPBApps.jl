
#---------------------------------------------------------------------
# this function performs the solution of the approximate factorization
# step in the y-direction for all five matrix components
# simultaneously. The Thomas algorithm is employed to solve the
# systems for the y-lines. Boundary conditions are non-periodic
#---------------------------------------------------------------------

function y_solve(_::Val{ncells}, # ::Int64,
                 successor, # ::Vector{Int64},
                 predecessor, # ::Vector{Int64},
                 slice, # ::Array{Int64,2},
                 cell_size, # ::Array{Int64,2},
                 cell_start, # ::Array{Int64,2},
                 cell_end, # ::Array{Int64,2},
                 cell_coord, # ::Array{Int64,2},
                 lhs, # ::OffsetArray{Float32, 5, Array{Float32, 5}},
                 rhs, # ::OffsetArray{Float32, 5, Array{Float32, 5}},
                 rho_i, # ::OffsetArray{Float32, 4, Array{Float32, 4}},
                 vs, # ::OffsetArray{Float32, 4, Array{Float32, 4}},
                 rhoq,
                 cv,
                speed, # ::OffsetArray{Float32, 4, Array{Float32, 4}},
                 con43, # ::Float32, 
                 c3c4, # ::Float32, 
                 c1c5, # ::Float32, 
                 c2dtty1,
                 dy3, # ::Float32, 
                 dy5, # ::Float32, 
                 dy1, # ::Float32, 
                 dymax, # ::Float32, 
                 dtty1,
                 dtty2, # ::Float32, 
                 comz5, # ::Float32, 
                 comz4, # ::Float32, 
                 comz1, # ::Float32, 
                 comz6, # ::Float32,
                 in_buffer, # ::Vector{Float32},
                 out_buffer, # ::Vector{Float32},
                 comm_solve, # ::MPI.Comm
                 requests,
                 s,
                 timeron
                 ) where ncells

#       requests = Array{MPI.Request}(undef,2)
#       s = Array{Float32}(undef,5)

#---------------------------------------------------------------------
# now do a sweep on a layer-by-layer basis, i.e. sweeping through cells
# on this node in the direction of increasing i for the forward sweep,
# and after that reversing the direction for the backsubstitution  
#---------------------------------------------------------------------

       if (timeron) timer_start(t_ysolve) end
#---------------------------------------------------------------------
#                          FORWARD ELIMINATION  
#---------------------------------------------------------------------
       for stage = 1:ncells
          c      = slice[2, stage]

          jstart = 0
          jend   = cell_size[2, c]-1

          isize     = cell_size[1, c]
          ksize     = cell_size[3, c]
          ip        = cell_coord[1, c]-1
          kp        = cell_coord[3, c]-1

          buffer_size = (isize-cell_start[1, c]-cell_end[1, c]) *(
                         ksize-cell_start[3, c]-cell_end[3, c])

          if stage != 1

#---------------------------------------------------------------------
#            if this is not the first processor in this row of cells, 
#            receive data from predecessor containing the right hand
#            sides and the upper diagonal elements of the previous two rows
#---------------------------------------------------------------------

             if (timeron) timer_start(t_ycomm) end
             requests[1] = MPI.Irecv!(in_buffer, predecessor[2], DEFAULT_TAG, comm_solve)
             if (timeron) timer_stop(t_ycomm) end

#---------------------------------------------------------------------
#            communication has already been started. 
#            compute the left hand side while waiting for the msg
#---------------------------------------------------------------------
             lhsy(c,
                  cell_size,
                  cell_start,
                  cell_end,
                  lhs,
                  rho_i,
                  vs,
                  rhoq,
                  cv,
                  speed,
                  con43, c3c4, c1c5, c2dtty1,
                  dy3, dy5, dy1, dymax, dtty1, dtty2, 
                  comz5, comz4, comz1, comz6)

#---------------------------------------------------------------------
#            wait for pending communication to complete
#            This waits on the current receive and on the send
#            from the previous stage. They always come in pairs. 
#---------------------------------------------------------------------
             if (timeron) timer_start(t_ycomm) end
             MPI.Waitall(requests)
             if (timeron) timer_stop(t_ycomm) end

#---------------------------------------------------------------------
#            unpack the buffer                                 
#---------------------------------------------------------------------
             j  = jstart
             j1 = jstart + 1
             n = 0
#---------------------------------------------------------------------
#            create a running pointer
#---------------------------------------------------------------------
             p = 0
             for k = cell_start[3, c]:ksize-cell_end[3, c]-1
                for i = cell_start[1, c]:isize-cell_end[1, c]-1
                   lhs[i, j, k, n+2, c] -= in_buffer[p+1] * lhs[i, j, k, n+1, c]
                   lhs[i, j, k, n+3, c] -= in_buffer[p+2] * lhs[i, j, k, n+1, c]
                   for m = 1:3
                      rhs[i, j, k, m, c] -= in_buffer[p+2+m] * lhs[i, j, k, n+1, c]
                   end
                   d            = in_buffer[p+6]
                   e            = in_buffer[p+7]
                   for m = 1:3
                      s[m] = in_buffer[p+7+m]
                   end
                   r1 = lhs[i, j, k, n+2, c]
                   lhs[i, j, k, n+3, c] -= d * r1
                   lhs[i, j, k, n+4, c] -= e * r1
                   for m = 1:3
                      rhs[i, j, k, m, c] -= s[m] * r1
                   end
                   r2 = lhs[i, j1, k, n+1, c]
                   lhs[i, j1, k, n+2, c] -= d * r2
                   lhs[i, j1, k, n+3, c] -= e * r2
                   for m = 1:3
                      rhs[i, j1, k, m, c] -= s[m] * r2
                   end
                   p = p + 10
                end
             end

             for m = 4:5
                n = (m-3)*5
                for k = cell_start[3, c]:ksize-cell_end[3, c]-1
                   for i = cell_start[1, c]:isize-cell_end[1, c]-1
                      lhs[i, j, k, n+2, c] -= in_buffer[p+1] * lhs[i, j, k, n+1, c]
                      lhs[i, j, k, n+3, c] -= in_buffer[p+2] * lhs[i, j, k, n+1, c]
                      rhs[i, j, k, m, c]   -= in_buffer[p+3] * lhs[i, j, k, n+1, c]
                      d                = in_buffer[p+4]
                      e                = in_buffer[p+5]
                      s[m]             = in_buffer[p+6]
                      r1 = lhs[i, j, k, n+2, c]
                      lhs[i, j, k, n+3, c] -= d * r1
                      lhs[i, j, k, n+4, c] -= e * r1
                      rhs[i, j, k, m, c]   -= s[m] * r1
                      r2 = lhs[i, j1, k, n+1, c]
                      lhs[i, j1, k, n+2, c] -= d * r2
                      lhs[i, j1, k, n+3, c] -= e * r2
                      rhs[i, j1, k, m, c]   -= s[m] * r2
                      p = p + 6
                   end
                end
             end

          else

#---------------------------------------------------------------------
#            if this IS the first cell, we still compute the lhs
#---------------------------------------------------------------------
             lhsy(c,
                  cell_size,
                  cell_start,
                  cell_end,
                  lhs,
                  rho_i,
                  vs,
                  rhoq,
                  cv,
                  speed,
                  con43, c3c4, c1c5, c2dtty1,
                  dy3, dy5, dy1, dymax, dtty1, dtty2, 
                  comz5, comz4, comz1, comz6)
               end

#---------------------------------------------------------------------
#         perform the Thomas algorithm; first, FORWARD ELIMINATION     
#---------------------------------------------------------------------
          n = 0

          for k = cell_start[3, c]:ksize-cell_end[3, c]-1
             for j = jstart:jend-2
                for i = cell_start[1, c]:isize-cell_end[1, c]-1
                   j1 = j  + 1
                   j2 = j  + 2
                   fac1 = 1.0e0/lhs[i, j, k, n+3, c]
                   lhs[i, j, k, n+4, c] *= fac1
                   lhs[i, j, k, n+5, c] *= fac1
                   for m = 1:3
                      rhs[i, j, k, m, c] *= fac1
                   end
                   lhs[i, j1, k, n+3, c] -= lhs[i, j1, k, n+2, c]*lhs[i, j, k, n+4, c]
                   lhs[i, j1, k, n+4, c] -= lhs[i, j1, k, n+2, c]*lhs[i, j, k, n+5, c]
                   for m = 1:3
                      rhs[i, j1, k, m, c] -= lhs[i, j1, k, n+2, c]*rhs[i, j, k, m, c]
                   end
                   lhs[i, j2, k, n+2, c] -= lhs[i, j2, k, n+1, c]*lhs[i, j, k, n+4, c]
                   lhs[i, j2, k, n+3, c] -= lhs[i, j2, k, n+1, c]*lhs[i, j, k, n+5, c]
                   for m = 1:3
                      rhs[i, j2, k, m, c] -= lhs[i, j2, k, n+1, c]*rhs[i, j, k, m, c]
                   end
                end
             end
          end

#---------------------------------------------------------------------
#         The last two rows in this grid block are a bit different, 
#         since they do not have two more rows available for the
#         elimination of off-diagonal entries
#---------------------------------------------------------------------

          j  = jend - 1
          j1 = jend
          for k = cell_start[3, c]:ksize-cell_end[3, c]-1
             for i = cell_start[1, c]:isize-cell_end[1, c]-1
                fac1 = 1.0e0/lhs[i, j, k, n+3, c]
                lhs[i, j, k, n+4, c] *= fac1
                lhs[i, j, k, n+5, c] *= fac1
                for m = 1:3
                   rhs[i, j, k, m, c] *= fac1
                end
                lhs[i, j1, k, n+3, c] -= lhs[i, j1, k, n+2, c]*lhs[i, j, k, n+4, c]
                lhs[i, j1, k, n+4, c] -= lhs[i, j1, k, n+2, c]*lhs[i, j, k, n+5, c]
                for m = 1:3
                   rhs[i, j1, k, m, c] -= lhs[i, j1, k, n+2, c]*rhs[i, j, k, m, c]
                end
#---------------------------------------------------------------------
#               scale the last row immediately (some of this is
#               overkill in case this is the last cell)
#---------------------------------------------------------------------
                fac2 = 1.0e0/lhs[i, j1, k, n+3, c]
                lhs[i, j1, k, n+4, c] *= fac2
                lhs[i, j1, k, n+5, c] *= fac2
                for m = 1:3
                   rhs[i, j1, k, m, c] *= fac2
                end
             end
          end

#---------------------------------------------------------------------
#         do the u+c and the u-c factors                 
#---------------------------------------------------------------------
          for m = 4:5
             n = (m-3)*5
             for k = cell_start[3, c]:ksize-cell_end[3, c]-1
                for j = jstart:jend-2
                   for i = cell_start[1, c]:isize-cell_end[1, c]-1
                   j1 = j  + 1
                   j2 = j  + 2
                   fac1 = 1.0e0/lhs[i, j, k, n+3, c]
                   lhs[i, j, k, n+4, c] *= fac1
                   lhs[i, j, k, n+5, c] *= fac1
                   rhs[i, j, k, m, c] *= fac1
                   lhs[i, j1, k, n+3, c] -= lhs[i, j1, k, n+2, c]*lhs[i, j, k, n+4, c]
                   lhs[i, j1, k, n+4, c] -= lhs[i, j1, k, n+2, c]*lhs[i, j, k, n+5, c]
                   rhs[i, j1, k, m, c] -= lhs[i, j1, k, n+2, c]*rhs[i, j, k, m, c]
                   lhs[i, j2, k, n+2, c]-=lhs[i, j2, k, n+1, c]*lhs[i, j, k, n+4, c]
                   lhs[i, j2, k, n+3, c] -= lhs[i, j2, k, n+1, c]*lhs[i, j, k, n+5, c]
                   rhs[i, j2, k, m, c] -= lhs[i, j2, k, n+1, c]*rhs[i, j, k, m, c]
                end
             end
          end

#---------------------------------------------------------------------
#            And again the last two rows separately
#---------------------------------------------------------------------
             j  = jend - 1
             j1 = jend
             for k = cell_start[3, c]:ksize-cell_end[3, c]-1
                for i = cell_start[1, c]:isize-cell_end[1, c]-1
                fac1 = 1.0e0/lhs[i, j, k, n+3, c]
                lhs[i, j, k, n+4, c] *= fac1
                lhs[i, j, k, n+5, c] *= fac1
                rhs[i, j, k, m, c] *= fac1
                lhs[i, j1, k, n+3, c] -= lhs[i, j1, k, n+2, c]*lhs[i, j, k, n+4, c]
                lhs[i, j1, k, n+4, c] -= lhs[i, j1, k, n+2, c]*lhs[i, j, k, n+5, c]
                rhs[i, j1, k, m, c]   -= lhs[i, j1, k, n+2, c]*rhs[i, j, k, m, c]
#---------------------------------------------------------------------
#               Scale the last row immediately (some of this is overkill
#               if this is the last cell)
#---------------------------------------------------------------------
                fac2 = 1.0e0/lhs[i, j1, k, n+3, c]
                lhs[i, j1, k, n+4, c] *= fac2
                lhs[i, j1, k, n+5, c] *= fac2
                rhs[i, j1, k, m, c]   *= fac2

             end
          end
       end

#---------------------------------------------------------------------
#         send information to the next processor, except when this
#         is the last grid block;
#---------------------------------------------------------------------

          if stage != ncells

#---------------------------------------------------------------------
#            create a running pointer for the send buffer  
#---------------------------------------------------------------------
             p = 0
             n = 0
             for k = cell_start[3, c]:ksize-cell_end[3, c]-1
                for i = cell_start[1, c]:isize-cell_end[1, c]-1
                   for j = jend-1:jend
                      out_buffer[p+1] = lhs[i, j, k, n+4, c]
                      out_buffer[p+2] = lhs[i, j, k, n+5, c]
                      for m = 1:3
                         out_buffer[p+2+m] = rhs[i, j, k, m, c]
                      end
                      p = p+5
                   end
                end
             end

             for m = 4:5
                n = (m-3)*5
                for k = cell_start[3, c]:ksize-cell_end[3, c]-1
                   for i = cell_start[1, c]:isize-cell_end[1, c]-1
                      for j = jend-1:jend
                         out_buffer[p+1] = lhs[i, j, k, n+4, c]
                         out_buffer[p+2] = lhs[i, j, k, n+5, c]
                         out_buffer[p+3] = rhs[i, j, k, m, c]
                         p = p + 3
                      end
                   end
                end
             end

#---------------------------------------------------------------------
#            pack and send the buffer
#---------------------------------------------------------------------
             if (timeron) timer_start(t_ycomm) end
             requests[2] = MPI.Isend(view(out_buffer,1:22*buffer_size), successor[2], DEFAULT_TAG, comm_solve)
             if (timeron) timer_stop(t_ycomm) end

          end
       end

#---------------------------------------------------------------------
#      now go in the reverse direction                      
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#                         BACKSUBSTITUTION 
#---------------------------------------------------------------------
       for stage = ncells:-1:1
          c = slice[2, stage]

          jstart = 0
          jend   = cell_size[2, c]-1

          isize = cell_size[1, c]
          ksize = cell_size[3, c]
          ip    = cell_coord[1, c]-1
          kp    = cell_coord[3, c]-1

          buffer_size = (isize-cell_start[1, c]-cell_end[1, c]) *(
                        ksize-cell_start[3, c]-cell_end[3, c])

          if stage != ncells

#---------------------------------------------------------------------
#            if this is not the starting cell in this row of cells, 
#            wait for a message to be received, containing the 
#            solution of the previous two stations     
#---------------------------------------------------------------------

             if (timeron) timer_start(t_ycomm) end
             requests[1] = MPI.Irecv!(in_buffer, successor[2], DEFAULT_TAG, comm_solve)
             if (timeron) timer_stop(t_ycomm) end


#---------------------------------------------------------------------
#            communication has already been started
#            while waiting, do the block-diagonal inversion for the 
#            cell that was just finished                
#---------------------------------------------------------------------

             pinvr(slice[2, stage+1],
                  cell_size,
                  cell_start,
                  cell_end,
                  rhs,
                  bt)

#---------------------------------------------------------------------
#            wait for pending communication to complete
#---------------------------------------------------------------------
             if (timeron) timer_start(t_ycomm) end
             #MPI.Waitall(2, requests, statuses, ERROR)
             MPI.Waitall(requests)
             if (timeron) timer_stop(t_ycomm) end

#---------------------------------------------------------------------
#            unpack the buffer for the first three factors         
#---------------------------------------------------------------------
             n = 0
             p = 0
             j  = jend
             j1 = j - 1
             for m = 1:3
                for k = cell_start[3, c]:ksize-cell_end[3, c]-1
                   for i = cell_start[1, c]:isize-cell_end[1, c]-1
                      sm1 = in_buffer[p+1]
                      sm2 = in_buffer[p+2]
                      rhs[i, j, k, m, c] -= lhs[i, j, k, n+4, c]*sm1 + lhs[i, j, k, n+5, c]*sm2
                      rhs[i, j1, k, m, c] -= lhs[i, j1, k, n+4, c] * rhs[i, j, k, m, c] + lhs[i, j1, k, n+5, c] * sm1
                      p = p + 2
                   end
                end
             end

#---------------------------------------------------------------------
#            now unpack the buffer for the remaining two factors
#---------------------------------------------------------------------
             for m = 4:5
                n = (m-3)*5
                for k = cell_start[3, c]:ksize-cell_end[3, c]-1
                   for i = cell_start[1, c]:isize-cell_end[1, c]-1
                      sm1 = in_buffer[p+1]
                      sm2 = in_buffer[p+2]
                      rhs[i, j, k, m, c] -= lhs[i, j, k, n+4, c]*sm1 + lhs[i, j, k, n+5, c]*sm2
                      rhs[i, j1, k, m, c] -= lhs[i, j1, k, n+4, c] * rhs[i, j, k, m, c] + lhs[i, j1, k, n+5, c] * sm1
                      p = p + 2
                   end
                end
             end

          else
#---------------------------------------------------------------------
#            now we know this is the first grid block on the back sweep,
#            so we don't need a message to start the substitution. 
#---------------------------------------------------------------------

             j  = jend - 1
             j1 = jend
             n = 0
             for m = 1:3
                for k = cell_start[3, c]:ksize-cell_end[3, c]-1
                   for i = cell_start[1, c]:isize-cell_end[1, c]-1
                      rhs[i, j, k, m, c] -= lhs[i, j, k, n+4, c]*rhs[i, j1, k, m, c]
                   end
                end
             end

             for m = 4:5
                n = (m-3)*5
                for k = cell_start[3, c]:ksize-cell_end[3, c]-1
                   for i = cell_start[1, c]:isize-cell_end[1, c]-1
                      rhs[i, j, k, m, c] -= lhs[i, j, k, n+4, c]*rhs[i, j1, k, m, c]
                   end
                end
             end
          end

#---------------------------------------------------------------------
#         Whether or not this is the last processor, we always have
#         to complete the back-substitution 
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#         The first three factors
#---------------------------------------------------------------------
          n = 0
          for m = 1:3
             for k = cell_start[3, c]:ksize-cell_end[3, c]-1
                for j = jend-2:-1:jstart
                   for i = cell_start[1, c]:isize-cell_end[1, c]-1
                      j1 = j + 1
                      j2 = j + 2
                      rhs[i, j, k, m, c] -= lhs[i, j, k, n+4, c]*rhs[i, j1, k, m, c] + lhs[i, j, k, n+5, c]*rhs[i, j2, k, m, c]
                   end
                end
             end
          end

#---------------------------------------------------------------------
#         And the remaining two
#---------------------------------------------------------------------
          for m = 4:5
             n = (m-3)*5
             for k = cell_start[3, c]:ksize-cell_end[3, c]-1
                for j = jend-2:-1:jstart
                   for i = cell_start[1, c]:isize-cell_end[1, c]-1
                      j1 = j  + 1
                      j2 = j1 + 1
                      rhs[i, j, k, m, c] -= lhs[i, j, k, n+4, c]*rhs[i, j1, k, m, c] + lhs[i, j, k, n+5, c]*rhs[i, j2, k, m, c]
                   end
                end
             end
          end

#---------------------------------------------------------------------
#         send on information to the previous processor, if needed
#---------------------------------------------------------------------
          if stage !=  1
             j  = jstart
             j1 = jstart + 1
             p = 0
             for m = 1:5
                for k = cell_start[3, c]:ksize-cell_end[3, c]-1
                   for i = cell_start[1, c]:isize-cell_end[1, c]-1
                      out_buffer[p+1] = rhs[i, j, k, m, c]
                      out_buffer[p+2] = rhs[i, j1, k, m, c]
                      p = p + 2
                   end
                end
             end

#---------------------------------------------------------------------
#            pack and send the buffer
#---------------------------------------------------------------------

             if (timeron) timer_start(t_ycomm) end
             requests[2] = MPI.Isend(view(out_buffer,1:10*buffer_size), predecessor[2], DEFAULT_TAG, comm_solve)
             if (timeron) timer_stop(t_ycomm) end

          end

#---------------------------------------------------------------------
#         If this was the last stage, do the block-diagonal inversion          
#---------------------------------------------------------------------
          if (stage == 1) 
            pinvr(c,
                  cell_size,
                  cell_start,
                  cell_end,
                  rhs,
                  bt) 
         end

       end

       if (timeron) timer_stop(t_ysolve) end

       return nothing
end