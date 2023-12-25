#---------------------------------------------------------------------
# this function performs the solution of the approximate factorization
# step in the z-direction for all five matrix components
# simultaneously. The Thomas algorithm is employed to solve the
# systems for the z-lines. Boundary conditions are non-periodic
#---------------------------------------------------------------------

function z_solve(_::Val{ncells}, # ::Int64,
                 successor, # ::Vector{Int64},
                 predecessor, # ::Vector{Int64},
                 slice, # ::Array{Int64,2},
                 cell_size, # ::Array{Int64,2},
                 cell_start, # ::Array{Int64,2},
                 cell_end, # ::Array{Int64,2},
                 cell_coord, # ::Array{Int64,2},
                 lhs, # ::OffsetArray{Float64, 5, Array{Float64, 5}},
                 rhs, # ::OffsetArray{Float64, 5, Array{Float64, 5}},
                 u,
                 us,
                 vs,
                 qs,
                 ainv,
                 rho_i, # ::OffsetArray{Float64, 4, Array{Float64, 4}},
                 ws, # ::OffsetArray{Float64, 4, Array{Float64, 4}},
                 rhos,
                 cv,
                 speed, # ::OffsetArray{Float64, 4, Array{Float64, 4}},
                 dz4, # ::Float64, 
                 dz5, # ::Float64,
                 dz1, # ::Float64, 
                 con43, # ::Float64, 
                 c3c4, # ::Float64, 
                 c1c5, # ::Float64, 
                 dttz2, # ::Float64, 
                 dttz1, # ::Float64, 
                 dzmax, # ::Float64,
                 comz5, # ::Float64, 
                 comz4, # ::Float64, 
                 comz1, # ::Float64, 
                 comz6, # ::Float64,
                 c2dttz1,
                 in_buffer, # ::Vector{Float64},
                 out_buffer, # ::Vector{Float64},
                 comm_solve, # ::MPI.Comm
                 requests,
                 s,
                 timeron
                 ) where ncells

#       requests = Array{MPI.Request}(undef,2)
#       s = Array{Float64}(undef,5)

#---------------------------------------------------------------------
# now do a sweep on a layer-by-layer basis, i.e. sweeping through cells
# on this node in the direction of increasing i for the forward sweep,
# and after that reversing the direction for the backsubstitution  
#---------------------------------------------------------------------

       if (timeron) timer_start(t_zsolve) end
#---------------------------------------------------------------------
#                          FORWARD ELIMINATION  
#---------------------------------------------------------------------
       for stage = 1:ncells
          c         = slice[3, stage]

          kstart = 0
          kend   = cell_size[3, c]-1

          isize     = cell_size[1, c]
          jsize     = cell_size[2, c]
          ip        = cell_coord[1, c]-1
          jp        = cell_coord[2, c]-1

          buffer_size = (isize-cell_start[1, c]-cell_end[1, c]) *(jsize-cell_start[2, c]-cell_end[2, c])

          if stage != 1

#---------------------------------------------------------------------
#            if this is not the first processor in this row of cells, 
#            receive data from predecessor containing the right hand
#            sides and the upper diagonal elements of the previous two rows
#---------------------------------------------------------------------

             if (timeron) timer_start(t_zcomm) end
             requests[1] = MPI.Irecv!(in_buffer, predecessor[3], DEFAULT_TAG, comm_solve)
             if (timeron) timer_stop(t_zcomm) end

#---------------------------------------------------------------------
#            communication has already been started. 
#            compute the left hand side while waiting for the msg
#---------------------------------------------------------------------
              lhsz(c,
                  cell_size,
                  cell_start,
                  cell_end,
                  lhs,
                  rho_i,
                  ws,
                  rhos,
                  cv,
                  speed,
                  dz4, dz5, dz1, con43, c3c4, c1c5, dttz2, dttz1, dzmax,
                  comz5, comz4, comz1, comz6, c2dttz1)

#---------------------------------------------------------------------
#            wait for pending communication to complete
#---------------------------------------------------------------------
             if (timeron) timer_start(t_zcomm) end
             MPI.Waitall(requests)
              if (timeron) timer_stop(t_zcomm) end

#---------------------------------------------------------------------
#            unpack the buffer                                 
#---------------------------------------------------------------------
             k  = kstart
             k1 = kstart + 1
             n = 0

#---------------------------------------------------------------------
#            create a running pointer
#---------------------------------------------------------------------
             p = 0
             for j = cell_start[2, c]:jsize-cell_end[2, c]-1
                for i = cell_start[1, c]:isize-cell_end[1, c]-1
                   lhs[i, j, k, n+2, c] = lhs[i, j, k, n+2, c] -
                             in_buffer[p+1] * lhs[i, j, k, n+1, c]
                   lhs[i, j, k, n+3, c] = lhs[i, j, k, n+3, c] -
                             in_buffer[p+2] * lhs[i, j, k, n+1, c]
                   for m = 1:3
                      rhs[i, j, k, m, c] = rhs[i, j, k, m, c] -
                             in_buffer[p+2+m] * lhs[i, j, k, n+1, c]
                   end
                   d            = in_buffer[p+6]
                   e            = in_buffer[p+7]
                   for m = 1:3
                      s[m] = in_buffer[p+7+m]
                   end
                   r1 = lhs[i, j, k, n+2, c]
                   lhs[i, j, k, n+3, c] = lhs[i, j, k, n+3, c] - d * r1
                   lhs[i, j, k, n+4, c] = lhs[i, j, k, n+4, c] - e * r1
                   for m = 1:3
                      rhs[i, j, k, m, c] = rhs[i, j, k, m, c] - s[m] * r1
                   end
                   r2 = lhs[i, j, k1, n+1, c]
                   lhs[i, j, k1, n+2, c] = lhs[i, j, k1, n+2, c] - d * r2
                   lhs[i, j, k1, n+3, c] = lhs[i, j, k1, n+3, c] - e * r2
                   for m = 1:3
                      rhs[i, j, k1, m, c] = rhs[i, j, k1, m, c] - s[m] * r2
                   end
                   p = p + 10
                end
             end

             for m = 4:5
                n = (m-3)*5
                for j = cell_start[2, c]:jsize-cell_end[2, c]-1
                   for i = cell_start[1, c]:isize-cell_end[1, c]-1
                      lhs[i, j, k, n+2, c] = lhs[i, j, k, n+2, c] -
                                in_buffer[p+1] * lhs[i, j, k, n+1, c]
                      lhs[i, j, k, n+3, c] = lhs[i, j, k, n+3, c] -
                                in_buffer[p+2] * lhs[i, j, k, n+1, c]
                      rhs[i, j, k, m, c]   = rhs[i, j, k, m, c] -
                                in_buffer[p+3] * lhs[i, j, k, n+1, c]
                      d                = in_buffer[p+4]
                      e                = in_buffer[p+5]
                      s[m]             = in_buffer[p+6]
                      r1 = lhs[i, j, k, n+2, c]
                      lhs[i, j, k, n+3, c] = lhs[i, j, k, n+3, c] - d * r1
                      lhs[i, j, k, n+4, c] = lhs[i, j, k, n+4, c] - e * r1
                      rhs[i, j, k, m, c]   = rhs[i, j, k, m, c] - s[m] * r1
                      r2 = lhs[i, j, k1, n+1, c]
                      lhs[i, j, k1, n+2, c] = lhs[i, j, k1, n+2, c] - d * r2
                      lhs[i, j, k1, n+3, c] = lhs[i, j, k1, n+3, c] - e * r2
                      rhs[i, j, k1, m, c]   = rhs[i, j, k1, m, c] - s[m] * r2
                      p = p + 6
                   end
                end
             end

          else

#---------------------------------------------------------------------
#            if this IS the first cell, we still compute the lhs
#---------------------------------------------------------------------
           lhsz(c,
                  cell_size,
                  cell_start,
                  cell_end,
                  lhs,
                  rho_i,
                  ws,
                  rhos,
                  cv,
                  speed,
                  dz4, dz5, dz1, con43, c3c4, c1c5, dttz2, dttz1, dzmax,
                  comz5, comz4, comz1, comz6, c2dttz1)
          end

#---------------------------------------------------------------------
#         perform the Thomas algorithm; first, FORWARD ELIMINATION     
#---------------------------------------------------------------------
          n = 0

          for k = kstart:kend-2
             for j = cell_start[2, c]:jsize-cell_end[2, c]-1
                for i = cell_start[1, c]:isize-cell_end[1, c]-1
                   k1 = k  + 1
                   k2 = k  + 2
                   fac1               = 1.0e0/lhs[i, j, k, n+3, c]
                   lhs[i, j, k, n+4, c]   = fac1*lhs[i, j, k, n+4, c]
                   lhs[i, j, k, n+5, c]   = fac1*lhs[i, j, k, n+5, c]
                   for m = 1:3
                      rhs[i, j, k, m, c] = fac1*rhs[i, j, k, m, c]
                   end
                   lhs[i, j, k1, n+3, c] = lhs[i, j, k1, n+3, c] -
                               lhs[i, j, k1, n+2, c]*lhs[i, j, k, n+4, c]
                   lhs[i, j, k1, n+4, c] = lhs[i, j, k1, n+4, c] -
                               lhs[i, j, k1, n+2, c]*lhs[i, j, k, n+5, c]
                   for m = 1:3
                      rhs[i, j, k1, m, c] = rhs[i, j, k1, m, c] -
                               lhs[i, j, k1, n+2, c]*rhs[i, j, k, m, c]
                   end
                   lhs[i, j, k2, n+2, c] = lhs[i, j, k2, n+2, c] -
                               lhs[i, j, k2, n+1, c]*lhs[i, j, k, n+4, c]
                   lhs[i, j, k2, n+3, c] = lhs[i, j, k2, n+3, c] -
                               lhs[i, j, k2, n+1, c]*lhs[i, j, k, n+5, c]
                   for m = 1:3
                      rhs[i, j, k2, m, c] = rhs[i, j, k2, m, c] -
                               lhs[i, j, k2, n+1, c]*rhs[i, j, k, m, c]
                   end
                end
             end
          end

#---------------------------------------------------------------------
#         The last two rows in this grid block are a bit different, 
#         since they do not have two more rows available for the
#         elimination of off-diagonal entries
#---------------------------------------------------------------------
          k  = kend - 1
          k1 = kend
          for j = cell_start[2, c]:jsize-cell_end[2, c]-1
             for i = cell_start[1, c]:isize-cell_end[1, c]-1
                fac1               = 1.0e0/lhs[i, j, k, n+3, c]
                lhs[i, j, k, n+4, c]   = fac1*lhs[i, j, k, n+4, c]
                lhs[i, j, k, n+5, c]   = fac1*lhs[i, j, k, n+5, c]
                for m = 1:3
                   rhs[i, j, k, m, c] = fac1*rhs[i, j, k, m, c]
                end
                lhs[i, j, k1, n+3, c] = lhs[i, j, k1, n+3, c] -
                            lhs[i, j, k1, n+2, c]*lhs[i, j, k, n+4, c]
                lhs[i, j, k1, n+4, c] = lhs[i, j, k1, n+4, c] -
                            lhs[i, j, k1, n+2, c]*lhs[i, j, k, n+5, c]
                for m = 1:3
                   rhs[i, j, k1, m, c] = rhs[i, j, k1, m, c] -
                            lhs[i, j, k1, n+2, c]*rhs[i, j, k, m, c]
                end
#---------------------------------------------------------------------
#               scale the last row immediately (some of this is
#               overkill in case this is the last cell)
#---------------------------------------------------------------------
                fac2               = 1.0e0/lhs[i, j, k1, n+3, c]
                lhs[i, j, k1, n+4, c] = fac2*lhs[i, j, k1, n+4, c]
                lhs[i, j, k1, n+5, c] = fac2*lhs[i, j, k1, n+5, c]
                for m = 1:3
                   rhs[i, j, k1, m, c] = fac2*rhs[i, j, k1, m, c]
                end
             end
          end

#---------------------------------------------------------------------
#         do the u+c and the u-c factors               
#---------------------------------------------------------------------
          for m = 4:5
             n = (m-3)*5
             for k = kstart:kend-2
                for j = cell_start[2, c]:jsize-cell_end[2, c]-1
                   for i = cell_start[1, c]:isize-cell_end[1, c]-1
                   k1 = k  + 1
                   k2 = k  + 2
                   fac1               = 1.0e0/lhs[i, j, k, n+3, c]
                   lhs[i, j, k, n+4, c]   = fac1*lhs[i, j, k, n+4, c]
                   lhs[i, j, k, n+5, c]   = fac1*lhs[i, j, k, n+5, c]
                   rhs[i, j, k, m, c] = fac1*rhs[i, j, k, m, c]
                   lhs[i, j, k1, n+3, c] = lhs[i, j, k1, n+3, c] -
                               lhs[i, j, k1, n+2, c]*lhs[i, j, k, n+4, c]
                   lhs[i, j, k1, n+4, c] = lhs[i, j, k1, n+4, c] -
                               lhs[i, j, k1, n+2, c]*lhs[i, j, k, n+5, c]
                   rhs[i, j, k1, m, c] = rhs[i, j, k1, m, c] -
                               lhs[i, j, k1, n+2, c]*rhs[i, j, k, m, c]
                   lhs[i, j, k2, n+2, c] = lhs[i, j, k2, n+2, c] -
                               lhs[i, j, k2, n+1, c]*lhs[i, j, k, n+4, c]
                   lhs[i, j, k2, n+3, c] = lhs[i, j, k2, n+3, c] -
                               lhs[i, j, k2, n+1, c]*lhs[i, j, k, n+5, c]
                   rhs[i, j, k2, m, c] = rhs[i, j, k2, m, c] -
                               lhs[i, j, k2, n+1, c]*rhs[i, j, k, m, c]
                end
             end
          end

#---------------------------------------------------------------------
#            And again the last two rows separately
#---------------------------------------------------------------------
             k  = kend - 1
             k1 = kend
             for j = cell_start[2, c]:jsize-cell_end[2, c]-1
                for i = cell_start[1, c]:isize-cell_end[1, c]-1
                fac1               = 1.0e0/lhs[i, j, k, n+3, c]
                lhs[i, j, k, n+4, c]   = fac1*lhs[i, j, k, n+4, c]
                lhs[i, j, k, n+5, c]   = fac1*lhs[i, j, k, n+5, c]
                rhs[i, j, k, m, c]     = fac1*rhs[i, j, k, m, c]
                lhs[i, j, k1, n+3, c] = lhs[i, j, k1, n+3, c] -
                            lhs[i, j, k1, n+2, c]*lhs[i, j, k, n+4, c]
                lhs[i, j, k1, n+4, c] = lhs[i, j, k1, n+4, c] -
                            lhs[i, j, k1, n+2, c]*lhs[i, j, k, n+5, c]
                rhs[i, j, k1, m, c]   = rhs[i, j, k1, m, c] -
                            lhs[i, j, k1, n+2, c]*rhs[i, j, k, m, c]
#---------------------------------------------------------------------
#               Scale the last row immediately (some of this is overkill
#               if this is the last cell)
#---------------------------------------------------------------------
                fac2               = 1.0e0/lhs[i, j, k1, n+3, c]
                lhs[i, j, k1, n+4, c] = fac2*lhs[i, j, k1, n+4, c]
                lhs[i, j, k1, n+5, c] = fac2*lhs[i, j, k1, n+5, c]
                rhs[i, j, k1, m, c]   = fac2*rhs[i, j, k1, m, c]

             end
          end
       end

#---------------------------------------------------------------------
#         send information to the next processor, except when this
#         is the last grid block,
#---------------------------------------------------------------------

          if stage != ncells

#---------------------------------------------------------------------
#            create a running pointer for the send buffer  
#---------------------------------------------------------------------
             p = 0
             n = 0
             for j = cell_start[2, c]:jsize-cell_end[2, c]-1
                for i = cell_start[1, c]:isize-cell_end[1, c]-1
                   for k = kend-1:kend
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
                for j = cell_start[2, c]:jsize-cell_end[2, c]-1
                   for i = cell_start[1, c]:isize-cell_end[1, c]-1
                      for k = kend-1:kend
                         out_buffer[p+1] = lhs[i, j, k, n+4, c]
                         out_buffer[p+2] = lhs[i, j, k, n+5, c]
                         out_buffer[p+3] = rhs[i, j, k, m, c]
                         p = p + 3
                      end
                   end
                end
             end


             if (timeron) timer_start(t_zcomm) end
             requests[2] = MPI.Isend(view(out_buffer,1:22*buffer_size), successor[3], DEFAULT_TAG, comm_solve)
             if (timeron) timer_stop(t_zcomm) end

          end
       end

#---------------------------------------------------------------------
#      now go in the reverse direction                      
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#                         BACKSUBSTITUTION 
#---------------------------------------------------------------------
       for stage = ncells:-1:1
          c = slice[3, stage]

          kstart = 0
          kend   = cell_size[3, c]-1

          isize     = cell_size[1, c]
          jsize     = cell_size[2, c]
          ip        = cell_coord[1, c]-1
          jp        = cell_coord[2, c]-1

          buffer_size = (isize-cell_start[1, c]-cell_end[1, c]) *(jsize-cell_start[2, c]-cell_end[2, c])

          if stage != ncells

#---------------------------------------------------------------------
#            if this is not the starting cell in this row of cells, 
#            wait for a message to be received, containing the 
#            solution of the previous two stations     
#---------------------------------------------------------------------

             if (timeron) timer_start(t_zcomm) end
             requests[1] = MPI.Irecv!(in_buffer, successor[3], DEFAULT_TAG, comm_solve)
             if (timeron) timer_stop(t_zcomm) end


#---------------------------------------------------------------------
#            communication has already been started
#            while waiting, do the  block-diagonal inversion for the 
#            cell that was just finished                
#---------------------------------------------------------------------

             tzetar(slice[3, stage+1],
                     cell_size,
                     cell_start,
                     cell_end,
                     rhs,
                     u,
                     us,
                     vs,
                     ws,
                     qs,
                     speed,
                     ainv,
                     bt,
                     c2iv)

#---------------------------------------------------------------------
#            wait for pending communication to complete
#---------------------------------------------------------------------
             if (timeron) timer_start(t_zcomm) end
             MPI.Waitall(requests)
             if (timeron) timer_stop(t_zcomm) end

#---------------------------------------------------------------------
#            unpack the buffer for the first three factors         
#---------------------------------------------------------------------
             n = 0
             p = 0
             k  = kend
             k1 = k - 1
             for m = 1:3
                for j = cell_start[2, c]:jsize-cell_end[2, c]-1
                   for i = cell_start[1, c]:isize-cell_end[1, c]-1
                      sm1 = in_buffer[p+1]
                      sm2 = in_buffer[p+2]
                      rhs[i, j, k, m, c] = rhs[i, j, k, m, c] -
                              lhs[i, j, k, n+4, c]*sm1 -
                              lhs[i, j, k, n+5, c]*sm2
                      rhs[i, j, k1, m, c] = rhs[i, j, k1, m, c] -
                              lhs[i, j, k1, n+4, c] * rhs[i, j, k, m, c] -
                              lhs[i, j, k1, n+5, c] * sm1
                      p = p + 2
                   end
                end
             end

#---------------------------------------------------------------------
#            now unpack the buffer for the remaining two factors
#---------------------------------------------------------------------
             for m = 4:5
                n = (m-3)*5
                for j = cell_start[2, c]:jsize-cell_end[2, c]-1
                   for i = cell_start[1, c]:isize-cell_end[1, c]-1
                      sm1 = in_buffer[p+1]
                      sm2 = in_buffer[p+2]
                      rhs[i, j, k, m, c] = rhs[i, j, k, m, c] -
                              lhs[i, j, k, n+4, c]*sm1 -
                              lhs[i, j, k, n+5, c]*sm2
                      rhs[i, j, k1, m, c] = rhs[i, j, k1, m, c] -
                              lhs[i, j, k1, n+4, c] * rhs[i, j, k, m, c] -
                              lhs[i, j, k1, n+5, c] * sm1
                      p = p + 2
                   end
                end
             end

          else

#---------------------------------------------------------------------
#            now we know this is the first grid block on the back sweep,
#            so we don't need a message to start the substitution. 
#---------------------------------------------------------------------

             k  = kend - 1
             k1 = kend
             n = 0
             for m = 1:3
                for j = cell_start[2, c]:jsize-cell_end[2, c]-1
                   for i = cell_start[1, c]:isize-cell_end[1, c]-1
                      rhs[i, j, k, m, c] = rhs[i, j, k, m, c] -
                                   lhs[i, j, k, n+4, c]*rhs[i, j, k1, m, c]
                   end
                end
             end

             for m = 4:5
                n = (m-3)*5
                for j = cell_start[2, c]:jsize-cell_end[2, c]-1
                   for i = cell_start[1, c]:isize-cell_end[1, c]-1
                      rhs[i, j, k, m, c] = rhs[i, j, k, m, c] -
                                   lhs[i, j, k, n+4, c]*rhs[i, j, k1, m, c]
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
             for k = kend-2:-1:kstart
                for j = cell_start[2, c]:jsize-cell_end[2, c]-1
                   for i = cell_start[1, c]:isize-cell_end[1, c]-1
                      k1 = k  + 1
                      k2 = k  + 2
                      rhs[i, j, k, m, c] = rhs[i, j, k, m, c] -
                                lhs[i, j, k, n+4, c]*rhs[i, j, k1, m, c] -
                                lhs[i, j, k, n+5, c]*rhs[i, j, k2, m, c]
                   end
                end
             end
          end

#---------------------------------------------------------------------
#         And the remaining two
#---------------------------------------------------------------------
          for m = 4:5
             n = (m-3)*5
             for k = kend-2:-1:kstart
                for j = cell_start[2, c]:jsize-cell_end[2, c]-1
                   for i = cell_start[1, c]:isize-cell_end[1, c]-1
                      k1 = k  + 1
                      k2 = k  + 2
                      rhs[i, j, k, m, c] = rhs[i, j, k, m, c] -
                                lhs[i, j, k, n+4, c]*rhs[i, j, k1, m, c] -
                                lhs[i, j, k, n+5, c]*rhs[i, j, k2, m, c]
                   end
                end
             end
          end

#---------------------------------------------------------------------
#         send on information to the previous processor, if needed
#---------------------------------------------------------------------
          if stage !=  1
             k  = kstart
             k1 = kstart + 1
             p = 0
             for m = 1:5
                for j = cell_start[2, c]:jsize-cell_end[2, c]-1
                   for i = cell_start[1, c]:isize-cell_end[1, c]-1
                      out_buffer[p+1] = rhs[i, j, k, m, c]
                      out_buffer[p+2] = rhs[i, j, k1, m, c]
                      p = p + 2
                   end
                end
             end

             if (timeron) timer_start(t_zcomm) end
             requests[2] = MPI.Isend(view(out_buffer,1:10*buffer_size), predecessor[3], DEFAULT_TAG, comm_solve)
             if (timeron) timer_stop(t_zcomm) end

          end

#---------------------------------------------------------------------
#         If this was the last stage, do the block-diagonal inversion
#---------------------------------------------------------------------
          if (stage == 1) 
            tzetar(c,
                  cell_size,
                  cell_start,
                  cell_end,
                  rhs,
                  u,
                  us,
                  vs,
                  ws,
                  qs,
                  speed,
                  ainv,
                  bt,
                  c2iv) 
         end

       end

       if (timeron) timer_stop(t_zsolve) end

       return nothing
end







