
#---------------------------------------------------------------------
#     Performs line solves in Y direction by first factoring
#     the block-tridiagonal matrix into an upper triangular matrix, 
#     and then performing back substitution to solve for the unknow
#     vectors of each line.  
#     
#     Make sure we treat elements zero to cell_size in the direction
#     of the sweep.
#---------------------------------------------------------------------

function y_solve(
     MAX_CELL_DIM,
     IMAX,
     JMAX,
     KMAX,
     cell_coord,
     cell_size,
     cell_start,
     cell_end,
     slice,
     u,
     rhs,
     lhsc,
     backsub_info,
     in_buffer,
     out_buffer,
     fjac,
     njac,
     lhsa,
     lhsb,
     qs,
     dt,
     timeron,
     ncells_v::Val{ncells},
     ty1,
     ty2,
     comm_solve,     
     predecessor,
     successor,
     utmp,
) where ncells

      jstart = 0

     if timeron timer_start(t_ysolve) end
#---------------------------------------------------------------------
#     in our terminology stage is the number of the cell in the y-direction
#     i.e. stage = 1 means the start of the line stage=ncells means end
#---------------------------------------------------------------------
      for stage = 1:ncells
         c = slice[2, stage]
         isize = cell_size[1, c] - 1
         jsize = cell_size[2, c] - 1
         ksize = cell_size[3, c] - 1

#---------------------------------------------------------------------
#     set last-cell flag
#---------------------------------------------------------------------
         if stage == ncells
            LAST = 1
         else
            LAST = 0
         end

         if stage == 1
#---------------------------------------------------------------------
#     This is the first cell, so solve without receiving data
#---------------------------------------------------------------------
            FIRST = 1
#            call lhsy(c)
            y_solve_cell(FIRST, LAST, c,
                         cell_size,
                         cell_start,
                         cell_end,    
                         u,
                         rhs,
                         lhsc,
                         fjac,
                         njac,
                         lhsa,
                         lhsb,
                         qs,
                         dt,
                         ty1,
                         ty2,
                         utmp,
                         )
         else
#---------------------------------------------------------------------
#     Not the first cell of this line, so receive info from
#     processor working on preceeding cell
#---------------------------------------------------------------------
            FIRST = 0
           if timeron timer_start(t_ycomm) end
            recv_id[] = y_receive_solve_info(c,
                                             MAX_CELL_DIM,
                                             cell_coord,
                                             out_buffer,
                                             ncells_v,
                                             comm_solve,
                                             predecessor
                                        )
#---------------------------------------------------------------------
#     overlap computations and communications
#---------------------------------------------------------------------
#            call lhsy(c)
#---------------------------------------------------------------------
#     wait for completion
#---------------------------------------------------------------------
            MPI.Wait(send_id[])
            MPI.Wait(recv_id[])
           if timeron timer_stop(t_ycomm) end
#---------------------------------------------------------------------
#     install C'(jstart+1) and rhs'(jstart+1) to be used in this cell
#---------------------------------------------------------------------
            y_unpack_solve_info(c,
                              IMAX,
                              KMAX,
                              rhs,
                              lhsc,
                              out_buffer,
                    )
            y_solve_cell(FIRST, LAST, c,
                         cell_size,
                         cell_start,
                         cell_end,    
                         u,
                         rhs,
                         lhsc,
                         fjac,
                         njac,
                         lhsa,
                         lhsb,
                         qs,
                         dt,
                         ty1,
                         ty2,
                         utmp,
                         )
         end

         if (LAST == 0) 
            send_id[] = y_send_solve_info(c,
                                        MAX_CELL_DIM,
                                        IMAX,
                                        KMAX,
                                        cell_coord,
                                        cell_size,
                                        rhs,
                                        lhsc,
                                        in_buffer,
                                        timeron,
                                        ncells_v,     
                                        comm_solve,    
                                        successor 
                                        ) 
         end
      end

#---------------------------------------------------------------------
#     now perform backsubstitution in reverse direction
#---------------------------------------------------------------------
      for stage = ncells:-1:1
         c = slice[2, stage]
         FIRST = 0
         LAST = 0
         if (stage == 1) FIRST = 1 end
         if stage == ncells
            LAST = 1
#---------------------------------------------------------------------
#     last cell, so perform back substitute without waiting
#---------------------------------------------------------------------
         y_backsubstitute(FIRST, LAST, c,
                              cell_size,
                              cell_start,
                              cell_end,
                              rhs,
                              lhsc,
                              backsub_info,
                              )
         else
           if timeron timer_start(t_ycomm) end
            recv_id[] = y_receive_backsub_info(c,
                                             MAX_CELL_DIM, 
                                             cell_coord,
                                             out_buffer,
                                             ncells_v,
                                             comm_solve,     
                                             successor
                                   )
            MPI.Wait(send_id[])
            MPI.Wait(recv_id[])
           if timeron timer_stop(t_ycomm) end
            y_unpack_backsub_info(c,
                                   IMAX,
                                   KMAX,
                                   backsub_info,
                                   out_buffer,
                         )
            y_backsubstitute(FIRST, LAST, c,
                              cell_size,
                              cell_start,
                              cell_end,
                              rhs,
                              lhsc,
                              backsub_info,
                              )
         end
         if (FIRST == 0) 
            send_id[] = y_send_backsub_info(c,
                                             MAX_CELL_DIM,
                                             IMAX,
                                             KMAX,
                                             cell_coord,
                                             rhs,
                                             in_buffer,
                                             timeron,
                                             ncells_v,
                                             comm_solve,
                                             predecessor
                                             ) 
         end
      end

     if timeron timer_stop(t_ysolve) end

      return nothing
end


#---------------------------------------------------------------------
#     unpack C'(-1) and rhs'(-1) for
#     all i and k
#---------------------------------------------------------------------

function y_unpack_solve_info(c,
                                   IMAX,
                                   KMAX,
                                   rhs,
                                   lhsc,
                                   out_buffer,
                                   )

      jstart = 0
      ptr = 0
      for k = 0:KMAX-1
         for i = 0:IMAX-1
            for m = 1:BLOCK_SIZE
               for n = 1:BLOCK_SIZE
                   lhsc[m, n, i, jstart-1, k, c] = out_buffer[ptr+n]
               end
               ptr = ptr+BLOCK_SIZE
            end
            for n = 1:BLOCK_SIZE
               rhs[n, i, jstart-1, k, c] = out_buffer[ptr+n]
            end
            ptr = ptr+BLOCK_SIZE
         end
      end

      return nothing
end


#---------------------------------------------------------------------
#     pack up and send C'(jend) and rhs'(jend) for
#     all i and k
#---------------------------------------------------------------------

function y_send_solve_info(c,
                                   MAX_CELL_DIM,
                                   IMAX,
                                   KMAX,
                                   cell_coord,
                                   cell_size,
                                   rhs,
                                   lhsc,
                                   in_buffer,
                                   timeron,
                                   ::Val{ncells},     
                                   comm_solve,     
                                   successor
                                   ) where ncells

      jsize = cell_size[2, c]-1
      ip = cell_coord[1, c] - 1
      kp = cell_coord[3, c] - 1
      buffer_size = MAX_CELL_DIM*MAX_CELL_DIM*(BLOCK_SIZE*BLOCK_SIZE + BLOCK_SIZE)

#---------------------------------------------------------------------
#     pack up buffer
#---------------------------------------------------------------------
      ptr = 0
      for k = 0:KMAX-1
         for i = 0:IMAX-1
            for m = 1:BLOCK_SIZE
               for n = 1:BLOCK_SIZE
                  in_buffer[ptr+n] =  lhsc[m, n, i, jsize, k, c]
               end
               ptr = ptr+BLOCK_SIZE
            end
            for n = 1:BLOCK_SIZE
               in_buffer[ptr+n] = rhs[n, i, jsize, k, c]
            end
            ptr = ptr+BLOCK_SIZE
         end
      end

#---------------------------------------------------------------------
#     send buffer 
#---------------------------------------------------------------------
     if timeron timer_start(t_ycomm) end
      send_id = MPI.Isend(view(in_buffer,1:buffer_size),successor[2],SOUTH+ip+kp*ncells, comm_solve)
     if timeron timer_stop(t_ycomm) end

      return send_id
end


#---------------------------------------------------------------------
#     pack up and send u[jstart] for all i and k
#---------------------------------------------------------------------

function y_send_backsub_info(c,
                                   MAX_CELL_DIM,
                                   IMAX,
                                   KMAX,
                                   cell_coord,
                                   rhs,
                                   in_buffer,
                                   timeron,
                                   ::Val{ncells},
                                   comm_solve,
                                   predecessor
                                   )  where ncells

#---------------------------------------------------------------------
#     Send element 0 to previous processor
#---------------------------------------------------------------------
      jstart = 0
      ip = cell_coord[1, c]-1
      kp = cell_coord[3, c]-1
      buffer_size = MAX_CELL_DIM*MAX_CELL_DIM*BLOCK_SIZE
      ptr = 0
      for k = 0:KMAX-1
         for i = 0:IMAX-1
            for n = 1:BLOCK_SIZE
               in_buffer[ptr+n] = rhs[n, i, jstart, k, c]
            end
            ptr = ptr+BLOCK_SIZE
         end
      end
     if timeron timer_start(t_ycomm) end
      send_id = MPI.Isend(view(in_buffer,1:buffer_size), predecessor[2], NORTH+ip+kp*ncells, comm_solve)
     if timeron timer_stop(t_ycomm) end

      return send_id
end


#---------------------------------------------------------------------
#     unpack u[jsize] for all i and k
#---------------------------------------------------------------------

function y_unpack_backsub_info(c,
                                        IMAX,
                                        KMAX,
                                        backsub_info,
                                        out_buffer,
     )

      ptr = 0
      for k = 0:KMAX-1
         for i = 0:IMAX-1
            for n = 1:BLOCK_SIZE
               backsub_info[n, i, k, c] = out_buffer[ptr+n]
            end
            ptr = ptr+BLOCK_SIZE
         end
      end

      return nothing
end


#---------------------------------------------------------------------
#     post mpi receives
#---------------------------------------------------------------------

function y_receive_backsub_info(c,
                                        MAX_CELL_DIM, 
                                        cell_coord,
                                        out_buffer,
                                        ::Val{ncells},
                                        comm_solve,     
                                        successor
          )  where ncells

      ip = cell_coord[1, c] - 1
      kp = cell_coord[3, c] - 1
      buffer_size = MAX_CELL_DIM*MAX_CELL_DIM*BLOCK_SIZE
      recv_id = MPI.Irecv!(view(out_buffer, 1:buffer_size), successor[2], NORTH+ip+kp*ncells, comm_solve)
      return recv_id
end


#---------------------------------------------------------------------
#     post mpi receives 
#---------------------------------------------------------------------

function y_receive_solve_info(c,
                                        MAX_CELL_DIM,
                                        cell_coord,
                                        out_buffer,
                                        ::Val{ncells},
                                        comm_solve,
                                        predecessor
                                    )  where ncells

      ip = cell_coord[1, c] - 1
      kp = cell_coord[3, c] - 1
      buffer_size = MAX_CELL_DIM*MAX_CELL_DIM*(BLOCK_SIZE*BLOCK_SIZE + BLOCK_SIZE)
      recv_id = MPI.Irecv!(view(out_buffer, 1:buffer_size), predecessor[2], SOUTH+ip+kp*ncells, comm_solve)

      return recv_id
end


#---------------------------------------------------------------------
#     back solve: if last cell, then generate u[jsize]=rhs[jsize]
#     else assume u[jsize] is loaded in un pack backsub_info
#     so just use it
#     after call u[jstart] will be sent to next cell
#---------------------------------------------------------------------

function y_backsubstitute(FIRST, LAST, c,
                                   cell_size,
                                   cell_start,
                                   cell_end,
                                   rhs,
                                   lhsc,
                                   backsub_info,
                                   )

      jstart = 0
      isize = cell_size[1, c]-cell_end[1, c]-1
      jsize = cell_size[2, c]-1
      ksize = cell_size[3, c]-cell_end[3, c]-1
      if LAST == 0
         for k = cell_start[3, c]:ksize
            for i = cell_start[1, c]:isize
#---------------------------------------------------------------------
#     u[jsize] uses info from previous cell if not last cell
#---------------------------------------------------------------------
               for m = 1:BLOCK_SIZE
                  for n = 1:BLOCK_SIZE
                     rhs[m, i, jsize, k, c] -= lhsc[m, n, i, jsize, k, c]*backsub_info[n, i, k, c]
                  end
               end
            end
         end
      end
      for k = cell_start[3, c]:ksize
         for j = jsize-1:-1:jstart
            for i = cell_start[1, c]:isize
               for m = 1:BLOCK_SIZE
                  for n = 1:BLOCK_SIZE
                     rhs[m, i, j, k, c] -= lhsc[m, n, i, j, k, c]*rhs[n, i, j+1, k, c]
                  end
               end
            end
         end
      end

      return nothing
end


#---------------------------------------------------------------------
#     performs guaussian elimination on this cell.
#     
#     assumes that unpacking routines for non-first cells 
#     preload C' and rhs' from previous cell.
#     
#     assumed send happens outside this routine, but that
#     c'(JMAX) and rhs'(JMAX) will be sent to next cell
#---------------------------------------------------------------------

function y_solve_cell(FIRST, LAST, c,
                    cell_size,
                    cell_start,
                    cell_end,    
                    u,
                    rhs,
                    lhsc,
                    fjac,
                    njac,
                    lhsa,
                    lhsb,
                    qs,
                    dt,
                    ty1,
                    ty2,
                    utmp,
          )

      
      
      jstart = 0
      isize = cell_size[1, c]-cell_end[1, c]-1
      jsize = cell_size[2, c]-1
      ksize = cell_size[3, c]-cell_end[3, c]-1

      lhsabinit(lhsa, lhsb, jsize)

      for k = cell_start[3, c]:ksize
         for i = cell_start[1, c]:isize

#---------------------------------------------------------------------
#     This function computes the left hand side for the three y-factors   
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     Compute the indices for storing the tri-diagonal matrix;
#     determine a (labeled f) and n jacobians for cell c
#---------------------------------------------------------------------
            for j = cell_start[2, c]-1:cell_size[2, c]-cell_end[2, c]
               utmp[1,j] = 1.0e0 / u[1, i, j, k, c]
               utmp[2,j] = u[2, i, j, k, c]
               utmp[3,j] = u[3, i, j, k, c]
               utmp[4,j] = u[4, i, j, k, c]
               utmp[5,j] = u[5, i, j, k, c]
               utmp[6,j] = qs[i, j, k, c]
            end

            for j = cell_start[2, c]-1:cell_size[2, c]-cell_end[2, c]

               tmp1 = utmp[1,j]
               tmp2 = tmp1 * tmp1
               tmp3 = tmp1 * tmp2

               fjac[1, 1, j] = 0.0e+00
               fjac[1, 2, j] = 0.0e+00
               fjac[1, 3, j] = 1.0e+00
               fjac[1, 4, j] = 0.0e+00
               fjac[1, 5, j] = 0.0e+00

               fjac[2, 1, j] = - ( utmp[2,j]*utmp[3,j] ) * tmp2
               fjac[2, 2, j] = utmp[3,j] * tmp1
               fjac[2, 3, j] = utmp[2,j] * tmp1
               fjac[2, 4, j] = 0.0e+00
               fjac[2, 5, j] = 0.0e+00

               fjac[3, 1, j] = - ( utmp[3,j]*utmp[3,j]*tmp2) + c2 * utmp[6,j]
               fjac[3, 2, j] = - c2 *  utmp[2,j] * tmp1
               fjac[3, 3, j] = ( 2.0e+00 - c2 ) * utmp[3,j] * tmp1
               fjac[3, 4, j] = - c2 * utmp[4,j] * tmp1
               fjac[3, 5, j] = c2

               fjac[4, 1, j] = - ( utmp[3,j]*utmp[4,j] ) * tmp2
               fjac[4, 2, j] = 0.0e+00
               fjac[4, 3, j] = utmp[4,j] * tmp1
               fjac[4, 4, j] = utmp[3,j] * tmp1
               fjac[4, 5, j] = 0.0e+00

               fjac[5, 1, j] = ( c2 * 2.0e0 * utmp[6,j] - c1 * utmp[5,j] * tmp1 ) * utmp[3,j] * tmp1
               fjac[5, 2, j] = - c2 * utmp[2,j]*utmp[3,j] * tmp2
               fjac[5, 3, j] = c1 * utmp[5,j] * tmp1 - c2 * ( utmp[6,j] + utmp[3,j]*utmp[3,j] * tmp2 )
               fjac[5, 4, j] = - c2 * ( utmp[3,j]*utmp[4,j] ) * tmp2
               fjac[5, 5, j] = c1 * utmp[3,j] * tmp1

               njac[1, 1, j] = 0.0e+00
               njac[1, 2, j] = 0.0e+00
               njac[1, 3, j] = 0.0e+00
               njac[1, 4, j] = 0.0e+00
               njac[1, 5, j] = 0.0e+00

               njac[2, 1, j] = - c3c4 * tmp2 * utmp[2,j]
               njac[2, 2, j] =   c3c4 * tmp1
               njac[2, 3, j] =   0.0e+00
               njac[2, 4, j] =   0.0e+00
               njac[2, 5, j] =   0.0e+00

               njac[3, 1, j] = - con43 * c3c4 * tmp2 * utmp[3,j]
               njac[3, 2, j] =   0.0e+00
               njac[3, 3, j] =   con43 * c3c4 * tmp1
               njac[3, 4, j] =   0.0e+00
               njac[3, 5, j] =   0.0e+00

               njac[4, 1, j] = - c3c4 * tmp2 * utmp[4,j]
               njac[4, 2, j] =   0.0e+00
               njac[4, 3, j] =   0.0e+00
               njac[4, 4, j] =   c3c4 * tmp1
               njac[4, 5, j] =   0.0e+00

               njac[5, 1, j] = - (  c3c4-
                     c1345 ) * tmp3 * (utmp[2,j]^2)-
                     ( con43 * c3c4-
                     c1345 ) * tmp3 * (utmp[3,j]^2)-
                     ( c3c4 - c1345 ) * tmp3 * (utmp[4,j]^2)-
                     c1345 * tmp2 * utmp[5,j]

               njac[5, 2, j] = (  c3c4 - c1345 ) * tmp2 * utmp[2,j]
               njac[5, 3, j] = ( con43 * c3c4 - c1345 ) * tmp2 * utmp[3,j]
               njac[5, 4, j] = ( c3c4 - c1345 ) * tmp2 * utmp[4,j]
               njac[5, 5, j] = ( c1345 ) * tmp1

            end

#---------------------------------------------------------------------
#     now joacobians set, so form left hand side in y direction
#---------------------------------------------------------------------
            for j = cell_start[2, c]:jsize-cell_end[2, c]

               tmp1 = dt * ty1
               tmp2 = dt * ty2

                lhsa[1, 1, j] = - tmp2 * fjac[1, 1, j-1] - tmp1 * njac[1, 1, j-1]- tmp1 * dy1
                lhsa[1, 2, j] = - tmp2 * fjac[1, 2, j-1] - tmp1 * njac[1, 2, j-1]
                lhsa[1, 3, j] = - tmp2 * fjac[1, 3, j-1] - tmp1 * njac[1, 3, j-1]
                lhsa[1, 4, j] = - tmp2 * fjac[1, 4, j-1] - tmp1 * njac[1, 4, j-1]
                lhsa[1, 5, j] = - tmp2 * fjac[1, 5, j-1] - tmp1 * njac[1, 5, j-1]

                lhsa[2, 1, j] = - tmp2 * fjac[2, 1, j-1] - tmp1 * njac[2, 1, j-1]
                lhsa[2, 2, j] = - tmp2 * fjac[2, 2, j-1] - tmp1 * njac[2, 2, j-1] - tmp1 * dy2
                lhsa[2, 3, j] = - tmp2 * fjac[2, 3, j-1] - tmp1 * njac[2, 3, j-1]
                lhsa[2, 4, j] = - tmp2 * fjac[2, 4, j-1] - tmp1 * njac[2, 4, j-1]
                lhsa[2, 5, j] = - tmp2 * fjac[2, 5, j-1] - tmp1 * njac[2, 5, j-1]

                lhsa[3, 1, j] = - tmp2 * fjac[3, 1, j-1] - tmp1 * njac[3, 1, j-1]
                lhsa[3, 2, j] = - tmp2 * fjac[3, 2, j-1] - tmp1 * njac[3, 2, j-1]
                lhsa[3, 3, j] = - tmp2 * fjac[3, 3, j-1] - tmp1 * njac[3, 3, j-1] - tmp1 * dy3
                lhsa[3, 4, j] = - tmp2 * fjac[3, 4, j-1] - tmp1 * njac[3, 4, j-1]
                lhsa[3, 5, j] = - tmp2 * fjac[3, 5, j-1] - tmp1 * njac[3, 5, j-1]

                lhsa[4, 1, j] = - tmp2 * fjac[4, 1, j-1] - tmp1 * njac[4, 1, j-1]
                lhsa[4, 2, j] = - tmp2 * fjac[4, 2, j-1] - tmp1 * njac[4, 2, j-1]
                lhsa[4, 3, j] = - tmp2 * fjac[4, 3, j-1] - tmp1 * njac[4, 3, j-1]
                lhsa[4, 4, j] = - tmp2 * fjac[4, 4, j-1] - tmp1 * njac[4, 4, j-1] - tmp1 * dy4
                lhsa[4, 5, j] = - tmp2 * fjac[4, 5, j-1] - tmp1 * njac[4, 5, j-1]

                lhsa[5, 1, j] = - tmp2 * fjac[5, 1, j-1] - tmp1 * njac[5, 1, j-1]
                lhsa[5, 2, j] = - tmp2 * fjac[5, 2, j-1] - tmp1 * njac[5, 2, j-1]
                lhsa[5, 3, j] = - tmp2 * fjac[5, 3, j-1] - tmp1 * njac[5, 3, j-1]
                lhsa[5, 4, j] = - tmp2 * fjac[5, 4, j-1] - tmp1 * njac[5, 4, j-1]
                lhsa[5, 5, j] = - tmp2 * fjac[5, 5, j-1] - tmp1 * njac[5, 5, j-1]- tmp1 * dy5

                lhsb[1, 1, j] = 1.0e+00 + tmp1 * 2.0e+00 * njac[1, 1, j] + tmp1 * 2.0e+00 * dy1
                lhsb[1, 2, j] = tmp1 * 2.0e+00 * njac[1, 2, j]
                lhsb[1, 3, j] = tmp1 * 2.0e+00 * njac[1, 3, j]
                lhsb[1, 4, j] = tmp1 * 2.0e+00 * njac[1, 4, j]
                lhsb[1, 5, j] = tmp1 * 2.0e+00 * njac[1, 5, j]

                lhsb[2, 1, j] = tmp1 * 2.0e+00 * njac[2, 1, j]
                lhsb[2, 2, j] = 1.0e+00 + tmp1 * 2.0e+00 * njac[2, 2, j] +  tmp1 * 2.0e+00 * dy2
                lhsb[2, 3, j] = tmp1 * 2.0e+00 * njac[2, 3, j]
                lhsb[2, 4, j] = tmp1 * 2.0e+00 * njac[2, 4, j]
                lhsb[2, 5, j] = tmp1 * 2.0e+00 * njac[2, 5, j]

                lhsb[3, 1, j] = tmp1 * 2.0e+00 * njac[3, 1, j]
                lhsb[3, 2, j] = tmp1 * 2.0e+00 * njac[3, 2, j]
                lhsb[3, 3, j] = 1.0e+00 + tmp1 * 2.0e+00 * njac[3, 3, j] + tmp1 * 2.0e+00 * dy3
                lhsb[3, 4, j] = tmp1 * 2.0e+00 * njac[3, 4, j]
                lhsb[3, 5, j] = tmp1 * 2.0e+00 * njac[3, 5, j]

                lhsb[4, 1, j] = tmp1 * 2.0e+00 * njac[4, 1, j]
                lhsb[4, 2, j] = tmp1 * 2.0e+00 * njac[4, 2, j]
                lhsb[4, 3, j] = tmp1 * 2.0e+00 * njac[4, 3, j]
                lhsb[4, 4, j] = 1.0e+00 + tmp1 * 2.0e+00 * njac[4, 4, j] + tmp1 * 2.0e+00 * dy4
                lhsb[4, 5, j] = tmp1 * 2.0e+00 * njac[4, 5, j]

                lhsb[5, 1, j] = tmp1 * 2.0e+00 * njac[5, 1, j]
                lhsb[5, 2, j] = tmp1 * 2.0e+00 * njac[5, 2, j]
                lhsb[5, 3, j] = tmp1 * 2.0e+00 * njac[5, 3, j]
                lhsb[5, 4, j] = tmp1 * 2.0e+00 * njac[5, 4, j]
                lhsb[5, 5, j] = 1.0e+00 + tmp1 * 2.0e+00 * njac[5, 5, j] + tmp1 * 2.0e+00 * dy5

                lhsc[1, 1, i, j, k, c] =  tmp2 * fjac[1, 1, j+1] -
                     tmp1 * njac[1, 1, j+1]-
                     tmp1 * dy1
                lhsc[1, 2, i, j, k, c] =  tmp2 * fjac[1, 2, j+1] -
                     tmp1 * njac[1, 2, j+1]
                lhsc[1, 3, i, j, k, c] =  tmp2 * fjac[1, 3, j+1] -
                     tmp1 * njac[1, 3, j+1]
                lhsc[1, 4, i, j, k, c] =  tmp2 * fjac[1, 4, j+1] -
                     tmp1 * njac[1, 4, j+1]
                lhsc[1, 5, i, j, k, c] =  tmp2 * fjac[1, 5, j+1] -
                     tmp1 * njac[1, 5, j+1]

                lhsc[2, 1, i, j, k, c] =  tmp2 * fjac[2, 1, j+1] -
                     tmp1 * njac[2, 1, j+1]
                lhsc[2, 2, i, j, k, c] =  tmp2 * fjac[2, 2, j+1] -
                     tmp1 * njac[2, 2, j+1]-
                     tmp1 * dy2
                lhsc[2, 3, i, j, k, c] =  tmp2 * fjac[2, 3, j+1] -
                     tmp1 * njac[2, 3, j+1]
                lhsc[2, 4, i, j, k, c] =  tmp2 * fjac[2, 4, j+1] -
                     tmp1 * njac[2, 4, j+1]
                lhsc[2, 5, i, j, k, c] =  tmp2 * fjac[2, 5, j+1] -
                     tmp1 * njac[2, 5, j+1]

                lhsc[3, 1, i, j, k, c] =  tmp2 * fjac[3, 1, j+1] -
                     tmp1 * njac[3, 1, j+1]
                lhsc[3, 2, i, j, k, c] =  tmp2 * fjac[3, 2, j+1] -
                     tmp1 * njac[3, 2, j+1]
                lhsc[3, 3, i, j, k, c] =  tmp2 * fjac[3, 3, j+1] -
                     tmp1 * njac[3, 3, j+1]-
                     tmp1 * dy3
                lhsc[3, 4, i, j, k, c] =  tmp2 * fjac[3, 4, j+1] -
                     tmp1 * njac[3, 4, j+1]
                lhsc[3, 5, i, j, k, c] =  tmp2 * fjac[3, 5, j+1] -
                     tmp1 * njac[3, 5, j+1]

                lhsc[4, 1, i, j, k, c] =  tmp2 * fjac[4, 1, j+1] -
                     tmp1 * njac[4, 1, j+1]
                lhsc[4, 2, i, j, k, c] =  tmp2 * fjac[4, 2, j+1] -
                     tmp1 * njac[4, 2, j+1]
                lhsc[4, 3, i, j, k, c] =  tmp2 * fjac[4, 3, j+1] -
                     tmp1 * njac[4, 3, j+1]
                lhsc[4, 4, i, j, k, c] =  tmp2 * fjac[4, 4, j+1] -
                     tmp1 * njac[4, 4, j+1]-
                     tmp1 * dy4
                lhsc[4, 5, i, j, k, c] =  tmp2 * fjac[4, 5, j+1] -
                     tmp1 * njac[4, 5, j+1]

                lhsc[5, 1, i, j, k, c] =  tmp2 * fjac[5, 1, j+1]-
                     tmp1 * njac[5, 1, j+1]
                lhsc[5, 2, i, j, k, c] =  tmp2 * fjac[5, 2, j+1]-
                     tmp1 * njac[5, 2, j+1]
                lhsc[5, 3, i, j, k, c] =  tmp2 * fjac[5, 3, j+1]-
                     tmp1 * njac[5, 3, j+1]
                lhsc[5, 4, i, j, k, c] =  tmp2 * fjac[5, 4, j+1]-
                     tmp1 * njac[5, 4, j+1]
                lhsc[5, 5, i, j, k, c] =  tmp2 * fjac[5, 5, j+1]-
                     tmp1 * njac[5, 5, j+1]-
                     tmp1 * dy5

            end


#---------------------------------------------------------------------
#     outer most do loops - sweeping in i direction
#---------------------------------------------------------------------
            if FIRST == 1

#---------------------------------------------------------------------
#     multiply c[i,jstart,k] by b_inverse and copy back to c
#     multiply rhs[jstart] by b_inverse(jstart) and copy to rhs
#---------------------------------------------------------------------
               binvcrhs(lhsb, jstart, lhsc, i, jstart, k, c, rhs, i, jstart, k, c)

            end

#---------------------------------------------------------------------
#     begin inner most do loop
#     do all the elements of the cell unless last 
#---------------------------------------------------------------------
            for j = jstart+FIRST:jsize-LAST

#---------------------------------------------------------------------
#     subtract A*lhs_vector(j-1) from lhs_vector(j)
#     
#     rhs[j] = rhs[j] - A*rhs[j-1]
#---------------------------------------------------------------------
#3               matvec_sub(view(lhsa, 1:5, 1:5, j), view(rhs, 1:5, i, j-1, k, c), view(rhs, 1:5, i, j, k, c))
               matvec_sub(lhsa, j, rhs, i, j-1, k, c, rhs, i, j, k, c)

#---------------------------------------------------------------------
#      b[j] =  b[j] - C(j-1)*A(j)
#---------------------------------------------------------------------
              matmul_sub2(view(lhsa, 1:5, 1:5, j), view(lhsc, 1:5, 1:5, i, j-1, k, c), view(lhsb, 1:5, 1:5, j))
#               matmul_sub(lhsa, j, lhsc, i, j-1, k, c, lhsb, j)

#---------------------------------------------------------------------
#     multiply c[i,j,k] by b_inverse and copy back to c
#     multiply rhs[i,1,k] by b_inverse(i,1,k) and copy to rhs
#---------------------------------------------------------------------
#               binvcrhs(view(lhsb, 1:5, 1:5, j), view(lhsc, 1:5, 1:5, i, j, k, c), view(rhs, 1:5, i, j, k, c) )
               binvcrhs(lhsb, j, lhsc, i, j, k, c, rhs, i, j, k, c)

            end

#---------------------------------------------------------------------
#     Now finish up special cases for last cell
#---------------------------------------------------------------------
            if LAST == 1

#---------------------------------------------------------------------
#     rhs[jsize] = rhs[jsize] - A*rhs[jsize-1]
#---------------------------------------------------------------------
#               matvec_sub(view(lhsa, 1:5, 1:5, jsize), view(rhs,1:5, i, jsize-1, k, c), view(rhs, 1:5, i, jsize, k, c))
               matvec_sub(lhsa, jsize, rhs, i, jsize-1, k, c, rhs, i, jsize, k, c)

#---------------------------------------------------------------------
#      b[jsize] =  b[jsize] - C(jsize-1)*A(jsize)
#     call matmul_sub(aa,i,jsize,k,c,
#     $              cc,i,jsize-1,k,c,bb,i,jsize,k,c)
#---------------------------------------------------------------------
               matmul_sub2(view(lhsa, 1:5, 1:5, jsize), view(lhsc, 1:5, 1:5, i, jsize-1, k, c), view(lhsb, 1:5, 1:5, jsize))
#               matmul_sub(lhsa, jsize, lhsc, i, jsize-1, k, c, lhsb, jsize)

#---------------------------------------------------------------------
#     multiply rhs[jsize] by b_inverse(jsize) and copy to rhs
#---------------------------------------------------------------------
#               binvrhs(view(lhsb,1:5, 1:5, jsize), view(rhs, 1:5, i, jsize, k, c))
               binvrhs(lhsb, jsize, rhs, i, jsize, k, c)

            end
         end
      end


      return nothing
end



