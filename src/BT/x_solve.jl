#---------------------------------------------------------------------
#     
#     Performs line solves in X direction by first factoring
#     the block-tridiagonal matrix into an upper triangular matrix, 
#     and then performing back substitution to solve for the unknow
#     vectors of each line.  
#     
#     Make sure we treat elements zero to cell_size in the direction
#     of the sweep.
#     
#---------------------------------------------------------------------


function x_solve(
               MAX_CELL_DIM,
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
               rho_i,
               dt,
               timeron,
               ncells_v::Val{ncells},
               tx1,
               tx2,
               comm_solve,
               predecessor,
               successor,
    ) where ncells
 @inbounds begin

      #istart = 0

     if timeron timer_start(t_xsolve) end
#---------------------------------------------------------------------
#     in our terminology stage is the number of the cell in the x-direction
#     i.e. stage = 1 means the start of the line stage=ncells means end
#---------------------------------------------------------------------
      for stage = 1:ncells
         c = slice[1, stage]
         #isize = cell_size[1, c] - 1
         #jsize = cell_size[2, c] - 1
         #ksize = cell_size[3, c] - 1

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
#            call lhsx(c)
            x_solve_cell(FIRST, LAST, c,
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
                         rho_i,
                         dt,
                         tx1,
                         tx2,
                         )
         else
#---------------------------------------------------------------------
#     Not the first cell of this line, so receive info from
#     processor working on preceeding cell
#---------------------------------------------------------------------
            FIRST = 0
           if timeron timer_start(t_xcomm) end
            recv_id[] = x_receive_solve_info(c,
                                             MAX_CELL_DIM,
                                             cell_coord,
                                             out_buffer,
                                             ncells_v,
                                             comm_solve,
                                             predecessor,
                                             )
#---------------------------------------------------------------------
#     overlap computations and communications
#---------------------------------------------------------------------
#            call lhsx(c)
#---------------------------------------------------------------------
#     wait for completion
#---------------------------------------------------------------------
            MPI.Wait(send_id[])
            MPI.Wait(recv_id[])

           if timeron timer_stop(t_xcomm) end
#---------------------------------------------------------------------
#     install C'(istart) and rhs'(istart) to be used in this cell
#---------------------------------------------------------------------
            x_unpack_solve_info(c,
                              JMAX,
                              KMAX,
                              rhs,
                              lhsc,
                              out_buffer,
                              )
            x_solve_cell(FIRST, LAST, c,
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
                         rho_i,
                         dt,
                         tx1,
                         tx2,
                         )
         end

         if (LAST == 0) 
          send_id[] = x_send_solve_info(c,
                                        MAX_CELL_DIM,
                                        JMAX,
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
         c = slice[1, stage]
         FIRST = 0
         LAST = 0
         if (stage == 1) FIRST = 1 end
         if stage == ncells
            LAST = 1
#---------------------------------------------------------------------
#     last cell, so perform back substitute without waiting
#---------------------------------------------------------------------
            x_backsubstitute(FIRST, LAST, c,
                              cell_size,
                              cell_start,
                              cell_end,
                              rhs,
                              lhsc,
                              backsub_info,     
                              )
         else
           if timeron timer_start(t_xcomm) end
            recv_id[] = x_receive_backsub_info(c,
                                             MAX_CELL_DIM,
                                             cell_coord,
                                             out_buffer,
                                             ncells_v,
                                             comm_solve,     
                                             successor
                                             )

            MPI.Wait(send_id[])
            MPI.Wait(recv_id[])
               
           if timeron timer_stop(t_xcomm) end
            x_unpack_backsub_info(c,
                                   JMAX,
                                   KMAX,
                                   backsub_info,
                                   out_buffer,
                                   )
            x_backsubstitute(FIRST, LAST, c,
                              cell_size,
                              cell_start,
                              cell_end,
                              rhs,
                              lhsc,
                              backsub_info,     
                              )
         end
         if (FIRST == 0) 
            send_id[] = x_send_backsub_info(c,
                                             MAX_CELL_DIM,
                                             JMAX,
                                             KMAX,
                                             cell_coord,
                                             rhs,
                                             in_buffer,
                                             timeron,
                                             ncells_v,
                                             comm_solve,
                                             predecessor,
                                        ) 
         end
      end

     if timeron timer_stop(t_xsolve) end

      return nothing
   end
end


#---------------------------------------------------------------------
#     unpack C'(-1) and rhs'(-1) for
#     all j and k
#---------------------------------------------------------------------

 function x_unpack_solve_info(c,
                              JMAX,
                              KMAX,
                              rhs,
                              lhsc,
                              out_buffer,
                             )
@inbounds begin
   
      istart = 0
      ptr = 0
      for k = 0:KMAX-1
         for j = 0:JMAX-1
            for m = 1:BLOCK_SIZE
               for n = 1:BLOCK_SIZE
                   lhsc[m, n, istart-1, j, k, c] = out_buffer[ptr+n]
               end
               ptr = ptr+BLOCK_SIZE
            end
            for n = 1:BLOCK_SIZE
               rhs[n, istart-1, j, k, c] = out_buffer[ptr+n]
            end
            ptr = ptr+BLOCK_SIZE
         end
      end

      return nothing
   end
end


#---------------------------------------------------------------------
#     pack up and send C'(iend) and rhs'(iend) for
#     all j and k
#---------------------------------------------------------------------

function x_send_solve_info(c,
                                   MAX_CELL_DIM,
                                   JMAX,
                                   KMAX,
                                   cell_coord,
                                   cell_size,
                                   rhs,
                                   lhsc,
                                   in_buffer,
                                   timeron,
                                   ::Val{ncells},
                                   comm_solve,    
                                   successor, 
                         ) where ncells
 @inbounds begin

      isize = cell_size[1, c]-1
      jp = cell_coord[2, c] - 1
      kp = cell_coord[3, c] - 1
      buffer_size = MAX_CELL_DIM*MAX_CELL_DIM*(BLOCK_SIZE*BLOCK_SIZE + BLOCK_SIZE)

#---------------------------------------------------------------------
#     pack up buffer
#---------------------------------------------------------------------
      ptr = 0
      for k = 0:KMAX-1
         for j = 0:JMAX-1
            for m = 1:BLOCK_SIZE
               for n = 1:BLOCK_SIZE
                  in_buffer[ptr+n] =  lhsc[m, n, isize, j, k, c]
               end
               ptr = ptr+BLOCK_SIZE
            end
            for n = 1:BLOCK_SIZE
               in_buffer[ptr+n] = rhs[n, isize, j, k, c]
            end
            ptr = ptr+BLOCK_SIZE
         end
      end

#---------------------------------------------------------------------
#     send buffer 
#---------------------------------------------------------------------
     if timeron timer_start(t_xcomm) end
      send_id = MPI.Isend(view(in_buffer,1:buffer_size), successor[1], WEST+jp+kp*ncells, comm_solve)
     if timeron timer_stop(t_xcomm) end

      return send_id
   end
end


#---------------------------------------------------------------------
#     pack up and send u[istart] for all j and k
#---------------------------------------------------------------------

function x_send_backsub_info(c,
                                   MAX_CELL_DIM,
                                   JMAX,
                                   KMAX,
                                   cell_coord,
                                   rhs,
                                   in_buffer,
                                   timeron,
                                   ::Val{ncells},
                                   comm_solve,
                                   predecessor,
                              ) where ncells
 @inbounds begin

#---------------------------------------------------------------------
#     Send element 0 to previous processor
#---------------------------------------------------------------------
      istart = 0
      jp = cell_coord[2, c]-1
      kp = cell_coord[3, c]-1
      buffer_size = MAX_CELL_DIM*MAX_CELL_DIM*BLOCK_SIZE
      ptr = 0
      for k = 0:KMAX-1
         for j = 0:JMAX-1
            for n = 1:BLOCK_SIZE
               in_buffer[ptr+n] = rhs[n, istart, j, k, c]
            end
            ptr = ptr+BLOCK_SIZE
         end
      end
     if timeron timer_start(t_xcomm) end
      send_id = MPI.Isend(view(in_buffer,1:buffer_size), predecessor[1], EAST+jp+kp*ncells, comm_solve)
     if timeron timer_stop(t_xcomm) end

      return send_id
   end
end


#---------------------------------------------------------------------
#     unpack u[isize] for all j and k
#---------------------------------------------------------------------

function x_unpack_backsub_info(c,
                                 JMAX,
                                 KMAX,
                                 backsub_info,
                                 out_buffer,
                                 )
 @inbounds begin

      ptr = 0
      for k = 0:KMAX-1
         for j = 0:JMAX-1
            for n = 1:BLOCK_SIZE
               backsub_info[n, j, k, c] = out_buffer[ptr+n]
            end
            ptr = ptr+BLOCK_SIZE
         end
      end

      return nothing
   end
end


#---------------------------------------------------------------------
#     post mpi receives
#---------------------------------------------------------------------

function x_receive_backsub_info(c,
                                 MAX_CELL_DIM,
                                 cell_coord,
                                 out_buffer,
                                 ::Val{ncells},
                                 comm_solve,     
                                 successor,
                                 ) where ncells
 @inbounds begin

      jp = cell_coord[2, c] - 1
      kp = cell_coord[3, c] - 1
      buffer_size = MAX_CELL_DIM*MAX_CELL_DIM*BLOCK_SIZE
      recv_id = MPI.Irecv!(view(out_buffer, 1:buffer_size), successor[1], EAST+jp+kp*ncells, comm_solve)

      return recv_id
   end
end


#---------------------------------------------------------------------
#     post mpi receives 
#---------------------------------------------------------------------

function x_receive_solve_info(c,
                              MAX_CELL_DIM,
                              cell_coord,
                              out_buffer,
                              ::Val{ncells},
                              comm_solve,
                              predecessor,
                              ) where ncells
 @inbounds begin

      jp = cell_coord[2, c] - 1
      kp = cell_coord[3, c] - 1
      buffer_size = MAX_CELL_DIM*MAX_CELL_DIM*(BLOCK_SIZE*BLOCK_SIZE + BLOCK_SIZE)
      recv_id = MPI.Irecv!(view(out_buffer, 1:buffer_size), predecessor[1], WEST+jp+kp*ncells, comm_solve)

      return recv_id
 end
end

#---------------------------------------------------------------------
#     back solve: if last cell, then generate u[isize]=rhs[isize]
#     else assume u[isize] is loaded in un pack backsub_info
#     so just use it
#     after call u[istart] will be sent to next cell
#---------------------------------------------------------------------

function x_backsubstitute(FIRST, LAST, c,
                           cell_size,
                           cell_start,
                           cell_end,
                           rhs,
                           lhsc,
                           backsub_info,     
                           )
 @inbounds begin

      istart = 0
      isize = cell_size[1, c]-1
      jsize = cell_size[2, c]-cell_end[2, c]-1
      ksize = cell_size[3, c]-cell_end[3, c]-1
      if LAST == 0
         for k = cell_start[3, c]:ksize
            for j = cell_start[2, c]:jsize
#---------------------------------------------------------------------
#     u[isize] uses info from previous cell if not last cell
#---------------------------------------------------------------------
               for m = 1:BLOCK_SIZE
                  for n = 1:BLOCK_SIZE
                     rhs[m, isize, j, k, c] -= lhsc[m, n, isize, j, k, c]*backsub_info[n, j, k, c]
#---------------------------------------------------------------------
#     rhs[m,isize,j,k,c] = rhs[m,isize,j,k,c] 
#     $                    -  lhsc[m,n,isize,j,k,c]*rhs[n,isize+1,j,k,c]
#---------------------------------------------------------------------
                  end
               end
            end
         end
      end
      for k = cell_start[3, c]:ksize
         for j = cell_start[2, c]:jsize
            for i = isize-1:-1:istart
               for m = 1:BLOCK_SIZE
                  for n = 1:BLOCK_SIZE
                     rhs[m, i, j, k, c] -= lhsc[m, n, i, j, k, c]*rhs[n, i+1, j, k, c]
                  end
               end
            end
         end
      end

      return nothing
   end
end


#---------------------------------------------------------------------
#     performs guaussian elimination on this cell.
#     
#     assumes that unpacking routines for non-first cells 
#     preload C' and rhs' from previous cell.
#     
#     assumed send happens outside this routine, but that
#     c'(IMAX) and rhs'(IMAX) will be sent to next cell
#---------------------------------------------------------------------

function x_solve_cell(FIRST, LAST, c,
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
                        rho_i,
                        dt,
                        tx1,
                        tx2,
                        )
   @inbounds begin

      istart = 0
      isize = cell_size[1, c]-1
      jsize = cell_size[2, c]-cell_end[2, c]-1
      ksize = cell_size[3, c]-cell_end[3, c]-1

      lhsabinit(lhsa, lhsb, isize)

      for k = cell_start[3, c]:ksize
         for j = cell_start[2, c]:jsize

#---------------------------------------------------------------------
#     This function computes the left hand side in the xi-direction
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     determine a (labeled f) and n jacobians for cell c
#---------------------------------------------------------------------
            for i = cell_start[1, c]-1:cell_size[1, c] - cell_end[1, c]

               tmp1 = rho_i[i, j, k, c]
               tmp2 = tmp1 * tmp1
               tmp3 = tmp1 * tmp2
#---------------------------------------------------------------------
#     
#---------------------------------------------------------------------
               fjac[1, 1, i] = 0.0e+00
               fjac[1, 2, i] = 1.0e+00
               fjac[1, 3, i] = 0.0e+00
               fjac[1, 4, i] = 0.0e+00
               fjac[1, 5, i] = 0.0e+00 

               fjac[2, 1, i] = -(u[2, i, j, k, c] * tmp2 * u[2, i, j, k, c]) + c2 * qs[i, j, k, c]
               fjac[2, 2, i] = ( 2.0e+00 - c2 ) * ( u[2, i, j, k, c] * tmp1 )
               fjac[2, 3, i] = - c2 * ( u[3, i, j, k, c] * tmp1 )
               fjac[2, 4, i] = - c2 * ( u[4, i, j, k, c] * tmp1 )
               fjac[2, 5, i] = c2

               fjac[3, 1, i] = - ( u[2, i, j, k, c]*u[3, i, j, k, c] ) * tmp2
               fjac[3, 2, i] = u[3, i, j, k, c] * tmp1
               fjac[3, 3, i] = u[2, i, j, k, c] * tmp1
               fjac[3, 4, i] = 0.0e+00
               fjac[3, 5, i] = 0.0e+00

               fjac[4, 1, i] = - ( u[2, i, j, k, c]*u[4, i, j, k, c] ) * tmp2
               fjac[4, 2, i] = u[4, i, j, k, c] * tmp1
               fjac[4, 3, i] = 0.0e+00
               fjac[4, 4, i] = u[2, i, j, k, c] * tmp1
               fjac[4, 5, i] = 0.0e+00

               fjac[5, 1, i] = ( c2 * 2.0e0 * qs[i, j, k, c] - c1 * ( u[5, i, j, k, c] * tmp1 ) ) * ( u[2, i, j, k, c] * tmp1 )
               fjac[5, 2, i] = c1 *  u[5, i, j, k, c] * tmp1 - c2 * ( u[2, i, j, k, c]*u[2, i, j, k, c] * tmp2 + qs[i, j, k, c] )
               fjac[5, 3, i] = - c2 * ( u[3, i, j, k, c]*u[2, i, j, k, c] ) * tmp2
               fjac[5, 4, i] = - c2 * ( u[4, i, j, k, c]*u[2, i, j, k, c] ) * tmp2
               fjac[5, 5, i] = c1 * ( u[2, i, j, k, c] * tmp1 )

               njac[1, 1, i] = 0.0e+00
               njac[1, 2, i] = 0.0e+00
               njac[1, 3, i] = 0.0e+00
               njac[1, 4, i] = 0.0e+00
               njac[1, 5, i] = 0.0e+00

               njac[2, 1, i] = - con43 * c3c4 * tmp2 * u[2, i, j, k, c]
               njac[2, 2, i] =   con43 * c3c4 * tmp1
               njac[2, 3, i] =   0.0e+00
               njac[2, 4, i] =   0.0e+00
               njac[2, 5, i] =   0.0e+00

               njac[3, 1, i] = - c3c4 * tmp2 * u[3, i, j, k, c]
               njac[3, 2, i] =   0.0e+00
               njac[3, 3, i] =   c3c4 * tmp1
               njac[3, 4, i] =   0.0e+00
               njac[3, 5, i] =   0.0e+00

               njac[4, 1, i] = - c3c4 * tmp2 * u[4, i, j, k, c]
               njac[4, 2, i] =   0.0e+00
               njac[4, 3, i] =   0.0e+00
               njac[4, 4, i] =   c3c4 * tmp1
               njac[4, 5, i] =   0.0e+00

               njac[5, 1, i] = - ( con43 * c3c4-
                     c1345 ) * tmp3 * (u[2, i, j, k, c]^2)-
                     ( c3c4 - c1345 ) * tmp3 * (u[3, i, j, k, c]^2)-
                     ( c3c4 - c1345 ) * tmp3 * (u[4, i, j, k, c]^2)-
                     c1345 * tmp2 * u[5, i, j, k, c]

               njac[5, 2, i] = ( con43 * c3c4 - c1345 ) * tmp2 * u[2, i, j, k, c]
               njac[5, 3, i] = ( c3c4 - c1345 ) * tmp2 * u[3, i, j, k, c]
               njac[5, 4, i] = ( c3c4 - c1345 ) * tmp2 * u[4, i, j, k, c]
               njac[5, 5, i] = ( c1345 ) * tmp1

            end
#---------------------------------------------------------------------
#     now jacobians set, so form left hand side in x direction
#---------------------------------------------------------------------
            for i = cell_start[1, c]:isize - cell_end[1, c]

               tmp1 = dt * tx1
               tmp2 = dt * tx2

                lhsa[1, 1, i] = - tmp2 * fjac[1, 1, i-1] - tmp1 * njac[1, 1, i-1] - tmp1 * dx1
                lhsa[1, 2, i] = - tmp2 * fjac[1, 2, i-1] - tmp1 * njac[1, 2, i-1]
                lhsa[1, 3, i] = - tmp2 * fjac[1, 3, i-1] - tmp1 * njac[1, 3, i-1]
                lhsa[1, 4, i] = - tmp2 * fjac[1, 4, i-1] - tmp1 * njac[1, 4, i-1]
                lhsa[1, 5, i] = - tmp2 * fjac[1, 5, i-1] - tmp1 * njac[1, 5, i-1]

                lhsa[2, 1, i] = - tmp2 * fjac[2, 1, i-1] - tmp1 * njac[2, 1, i-1]
                lhsa[2, 2, i] = - tmp2 * fjac[2, 2, i-1] - tmp1 * njac[2, 2, i-1] - tmp1 * dx2
                lhsa[2, 3, i] = - tmp2 * fjac[2, 3, i-1] - tmp1 * njac[2, 3, i-1]
                lhsa[2, 4, i] = - tmp2 * fjac[2, 4, i-1] - tmp1 * njac[2, 4, i-1]
                lhsa[2, 5, i] = - tmp2 * fjac[2, 5, i-1] - tmp1 * njac[2, 5, i-1]

                lhsa[3, 1, i] = - tmp2 * fjac[3, 1, i-1] - tmp1 * njac[3, 1, i-1]
                lhsa[3, 2, i] = - tmp2 * fjac[3, 2, i-1] - tmp1 * njac[3, 2, i-1]
                lhsa[3, 3, i] = - tmp2 * fjac[3, 3, i-1] - tmp1 * njac[3, 3, i-1] - tmp1 * dx3
                lhsa[3, 4, i] = - tmp2 * fjac[3, 4, i-1] - tmp1 * njac[3, 4, i-1]
                lhsa[3, 5, i] = - tmp2 * fjac[3, 5, i-1] - tmp1 * njac[3, 5, i-1]

                lhsa[4, 1, i] = - tmp2 * fjac[4, 1, i-1] - tmp1 * njac[4, 1, i-1]
                lhsa[4, 2, i] = - tmp2 * fjac[4, 2, i-1] - tmp1 * njac[4, 2, i-1]
                lhsa[4, 3, i] = - tmp2 * fjac[4, 3, i-1] - tmp1 * njac[4, 3, i-1]
                lhsa[4, 4, i] = - tmp2 * fjac[4, 4, i-1] - tmp1 * njac[4, 4, i-1] - tmp1 * dx4
                lhsa[4, 5, i] = - tmp2 * fjac[4, 5, i-1] - tmp1 * njac[4, 5, i-1]

                lhsa[5, 1, i] = - tmp2 * fjac[5, 1, i-1] - tmp1 * njac[5, 1, i-1]
                lhsa[5, 2, i] = - tmp2 * fjac[5, 2, i-1] - tmp1 * njac[5, 2, i-1]
                lhsa[5, 3, i] = - tmp2 * fjac[5, 3, i-1] - tmp1 * njac[5, 3, i-1]
                lhsa[5, 4, i] = - tmp2 * fjac[5, 4, i-1] - tmp1 * njac[5, 4, i-1]
                lhsa[5, 5, i] = - tmp2 * fjac[5, 5, i-1] - tmp1 * njac[5, 5, i-1] - tmp1 * dx5

                lhsb[1, 1, i] = 1.0e+00 + tmp1 * 2.0e+00 * njac[1, 1, i] + tmp1 * 2.0e+00 * dx1
                lhsb[1, 2, i] = tmp1 * 2.0e+00 * njac[1, 2, i]
                lhsb[1, 3, i] = tmp1 * 2.0e+00 * njac[1, 3, i]
                lhsb[1, 4, i] = tmp1 * 2.0e+00 * njac[1, 4, i]
                lhsb[1, 5, i] = tmp1 * 2.0e+00 * njac[1, 5, i]

                lhsb[2, 1, i] = tmp1 * 2.0e+00 * njac[2, 1, i]
                lhsb[2, 2, i] = 1.0e+00 + tmp1 * 2.0e+00 * njac[2, 2, i] + tmp1 * 2.0e+00 * dx2
                lhsb[2, 3, i] = tmp1 * 2.0e+00 * njac[2, 3, i]
                lhsb[2, 4, i] = tmp1 * 2.0e+00 * njac[2, 4, i]
                lhsb[2, 5, i] = tmp1 * 2.0e+00 * njac[2, 5, i]

                lhsb[3, 1, i] = tmp1 * 2.0e+00 * njac[3, 1, i]
                lhsb[3, 2, i] = tmp1 * 2.0e+00 * njac[3, 2, i]
                lhsb[3, 3, i] = 1.0e+00 + tmp1 * 2.0e+00 * njac[3, 3, i] +tmp1 * 2.0e+00 * dx3
                lhsb[3, 4, i] = tmp1 * 2.0e+00 * njac[3, 4, i]
                lhsb[3, 5, i] = tmp1 * 2.0e+00 * njac[3, 5, i]

                lhsb[4, 1, i] = tmp1 * 2.0e+00 * njac[4, 1, i]
                lhsb[4, 2, i] = tmp1 * 2.0e+00 * njac[4, 2, i]
                lhsb[4, 3, i] = tmp1 * 2.0e+00 * njac[4, 3, i]
                lhsb[4, 4, i] = 1.0e+00 + tmp1 * 2.0e+00 * njac[4, 4, i] + tmp1 * 2.0e+00 * dx4
                lhsb[4, 5, i] = tmp1 * 2.0e+00 * njac[4, 5, i]

                lhsb[5, 1, i] = tmp1 * 2.0e+00 * njac[5, 1, i]
                lhsb[5, 2, i] = tmp1 * 2.0e+00 * njac[5, 2, i]
                lhsb[5, 3, i] = tmp1 * 2.0e+00 * njac[5, 3, i]
                lhsb[5, 4, i] = tmp1 * 2.0e+00 * njac[5, 4, i]
                lhsb[5, 5, i] = 1.0e+00 + tmp1 * 2.0e+00 * njac[5, 5, i] + tmp1 * 2.0e+00 * dx5

                lhsc[1, 1, i, j, k, c] =  tmp2 * fjac[1, 1, i+1] - tmp1 * njac[1, 1, i+1] - tmp1 * dx1
                lhsc[1, 2, i, j, k, c] =  tmp2 * fjac[1, 2, i+1] - tmp1 * njac[1, 2, i+1]
                lhsc[1, 3, i, j, k, c] =  tmp2 * fjac[1, 3, i+1] - tmp1 * njac[1, 3, i+1]
                lhsc[1, 4, i, j, k, c] =  tmp2 * fjac[1, 4, i+1] - tmp1 * njac[1, 4, i+1]
                lhsc[1, 5, i, j, k, c] =  tmp2 * fjac[1, 5, i+1] - tmp1 * njac[1, 5, i+1]

                lhsc[2, 1, i, j, k, c] =  tmp2 * fjac[2, 1, i+1] - tmp1 * njac[2, 1, i+1]
                lhsc[2, 2, i, j, k, c] =  tmp2 * fjac[2, 2, i+1] - tmp1 * njac[2, 2, i+1] - tmp1 * dx2
                lhsc[2, 3, i, j, k, c] =  tmp2 * fjac[2, 3, i+1] - tmp1 * njac[2, 3, i+1]
                lhsc[2, 4, i, j, k, c] =  tmp2 * fjac[2, 4, i+1] - tmp1 * njac[2, 4, i+1]
                lhsc[2, 5, i, j, k, c] =  tmp2 * fjac[2, 5, i+1] - tmp1 * njac[2, 5, i+1]

                lhsc[3, 1, i, j, k, c] =  tmp2 * fjac[3, 1, i+1] - tmp1 * njac[3, 1, i+1]
                lhsc[3, 2, i, j, k, c] =  tmp2 * fjac[3, 2, i+1] - tmp1 * njac[3, 2, i+1]
                lhsc[3, 3, i, j, k, c] =  tmp2 * fjac[3, 3, i+1] - tmp1 * njac[3, 3, i+1] - tmp1 * dx3
                lhsc[3, 4, i, j, k, c] =  tmp2 * fjac[3, 4, i+1] - tmp1 * njac[3, 4, i+1]
                lhsc[3, 5, i, j, k, c] =  tmp2 * fjac[3, 5, i+1] - tmp1 * njac[3, 5, i+1]

                lhsc[4, 1, i, j, k, c] =  tmp2 * fjac[4, 1, i+1] - tmp1 * njac[4, 1, i+1]
                lhsc[4, 2, i, j, k, c] =  tmp2 * fjac[4, 2, i+1] - tmp1 * njac[4, 2, i+1]
                lhsc[4, 3, i, j, k, c] =  tmp2 * fjac[4, 3, i+1] - tmp1 * njac[4, 3, i+1]
                lhsc[4, 4, i, j, k, c] =  tmp2 * fjac[4, 4, i+1] - tmp1 * njac[4, 4, i+1] - tmp1 * dx4
                lhsc[4, 5, i, j, k, c] =  tmp2 * fjac[4, 5, i+1] - tmp1 * njac[4, 5, i+1]

                lhsc[5, 1, i, j, k, c] =  tmp2 * fjac[5, 1, i+1] - tmp1 * njac[5, 1, i+1]
                lhsc[5, 2, i, j, k, c] =  tmp2 * fjac[5, 2, i+1] - tmp1 * njac[5, 2, i+1]
                lhsc[5, 3, i, j, k, c] =  tmp2 * fjac[5, 3, i+1] - tmp1 * njac[5, 3, i+1]
                lhsc[5, 4, i, j, k, c] =  tmp2 * fjac[5, 4, i+1] - tmp1 * njac[5, 4, i+1]
                lhsc[5, 5, i, j, k, c] =  tmp2 * fjac[5, 5, i+1] - tmp1 * njac[5, 5, i+1] - tmp1 * dx5

            end


#---------------------------------------------------------------------
#     outer most do loops - sweeping in i direction
#---------------------------------------------------------------------
            if FIRST == 1

#---------------------------------------------------------------------
#     multiply c[istart,j,k] by b_inverse and copy back to c
#     multiply rhs[istart] by b_inverse(istart) and copy to rhs
#---------------------------------------------------------------------
               binvcrhs(lhsb, istart, lhsc, istart, j, k, c, rhs, istart, j, k, c)

              #= i = istart
               ic1, ic2, ic3, ic4 = istart, j, k, c
               ir1, ir2, ir3, ir4 = istart, j, k, c

               pivot = 1.00e0/lhsb[1,1,i]
               lhsb[1,2,i] *= pivot
               lhsb[1,3,i] *= pivot
               lhsb[1,4,i] *= pivot
               lhsb[1,5,i] *= pivot
               lhsc[1,1,ic1,ic2,ic3,ic4] *= pivot
               lhsc[1,2,ic1,ic2,ic3,ic4] *= pivot
               lhsc[1,3,ic1,ic2,ic3,ic4] *= pivot
               lhsc[1,4,ic1,ic2,ic3,ic4] *= pivot
               lhsc[1,5,ic1,ic2,ic3,ic4] *= pivot
               rhs[1,ir1,ir2,ir3,ir4]   *= pivot
         
               coeff = lhsb[2,1,i]
               lhsb[2,2,i] -= coeff*lhsb[1,2,i]
               lhsb[2,3,i] -= coeff*lhsb[1,3,i]
               lhsb[2,4,i] -= coeff*lhsb[1,4,i]
               lhsb[2,5,i] -= coeff*lhsb[1,5,i]
               lhsc[2,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,1,ic1,ic2,ic3,ic4]
               lhsc[2,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,2,ic1,ic2,ic3,ic4]
               lhsc[2,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,3,ic1,ic2,ic3,ic4]
               lhsc[2,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,4,ic1,ic2,ic3,ic4]
               lhsc[2,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,5,ic1,ic2,ic3,ic4]
               rhs[2,ir1,ir2,ir3,ir4] -= coeff*rhs[1,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[3,1,i]
               lhsb[3,2,i] -= coeff*lhsb[1,2,i]
               lhsb[3,3,i] -= coeff*lhsb[1,3,i]
               lhsb[3,4,i] -= coeff*lhsb[1,4,i]
               lhsb[3,5,i] -= coeff*lhsb[1,5,i]
               lhsc[3,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,1,ic1,ic2,ic3,ic4]
               lhsc[3,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,2,ic1,ic2,ic3,ic4]
               lhsc[3,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,3,ic1,ic2,ic3,ic4]
               lhsc[3,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,4,ic1,ic2,ic3,ic4]
               lhsc[3,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,5,ic1,ic2,ic3,ic4]
               rhs[3,ir1,ir2,ir3,ir4] -= coeff*rhs[1,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[4,1,i]
               lhsb[4,2,i] -= coeff*lhsb[1,2,i]
               lhsb[4,3,i] -= coeff*lhsb[1,3,i]
               lhsb[4,4,i] -= coeff*lhsb[1,4,i]
               lhsb[4,5,i] -= coeff*lhsb[1,5,i]
               lhsc[4,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,1,ic1,ic2,ic3,ic4]
               lhsc[4,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,2,ic1,ic2,ic3,ic4]
               lhsc[4,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,3,ic1,ic2,ic3,ic4]
               lhsc[4,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,4,ic1,ic2,ic3,ic4]
               lhsc[4,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,5,ic1,ic2,ic3,ic4]
               rhs[4,ir1,ir2,ir3,ir4] -= coeff*rhs[1,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[5,1,i]
               lhsb[5,2,i] -= coeff*lhsb[1,2,i]
               lhsb[5,3,i] -= coeff*lhsb[1,3,i]
               lhsb[5,4,i] -= coeff*lhsb[1,4,i]
               lhsb[5,5,i] -= coeff*lhsb[1,5,i]
               lhsc[5,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,1,ic1,ic2,ic3,ic4]
               lhsc[5,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,2,ic1,ic2,ic3,ic4]
               lhsc[5,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,3,ic1,ic2,ic3,ic4]
               lhsc[5,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,4,ic1,ic2,ic3,ic4]
               lhsc[5,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,5,ic1,ic2,ic3,ic4]
               rhs[5,ir1,ir2,ir3,ir4] -= coeff*rhs[1,ir1,ir2,ir3,ir4]
         
         
               pivot = 1.00e0/lhsb[2,2,i]
               lhsb[2,3,i] *= pivot
               lhsb[2,4,i] *= pivot
               lhsb[2,5,i] *= pivot
               lhsc[2,1,ic1,ic2,ic3,ic4] *= pivot
               lhsc[2,2,ic1,ic2,ic3,ic4] *= pivot
               lhsc[2,3,ic1,ic2,ic3,ic4] *= pivot
               lhsc[2,4,ic1,ic2,ic3,ic4] *= pivot
               lhsc[2,5,ic1,ic2,ic3,ic4] *= pivot
               rhs[2,ir1,ir2,ir3,ir4] *= pivot
         
               coeff = lhsb[1,2,i]
               lhsb[1,3,i] -= coeff*lhsb[2,3,i]
               lhsb[1,4,i] -= coeff*lhsb[2,4,i]
               lhsb[1,5,i] -= coeff*lhsb[2,5,i]
               lhsc[1,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,1,ic1,ic2,ic3,ic4]
               lhsc[1,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,2,ic1,ic2,ic3,ic4]
               lhsc[1,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,3,ic1,ic2,ic3,ic4]
               lhsc[1,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,4,ic1,ic2,ic3,ic4]
               lhsc[1,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,5,ic1,ic2,ic3,ic4]
               rhs[1,ir1,ir2,ir3,ir4] -= coeff*rhs[2,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[3,2,i]
               lhsb[3,3,i] -= coeff*lhsb[2,3,i]
               lhsb[3,4,i] -= coeff*lhsb[2,4,i]
               lhsb[3,5,i] -= coeff*lhsb[2,5,i]
               lhsc[3,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,1,ic1,ic2,ic3,ic4]
               lhsc[3,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,2,ic1,ic2,ic3,ic4]
               lhsc[3,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,3,ic1,ic2,ic3,ic4]
               lhsc[3,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,4,ic1,ic2,ic3,ic4]
               lhsc[3,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,5,ic1,ic2,ic3,ic4]
               rhs[3,ir1,ir2,ir3,ir4] -= coeff*rhs[2,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[4,2,i]
               lhsb[4,3,i] -= coeff*lhsb[2,3,i]
               lhsb[4,4,i] -= coeff*lhsb[2,4,i]
               lhsb[4,5,i] -= coeff*lhsb[2,5,i]
               lhsc[4,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,1,ic1,ic2,ic3,ic4]
               lhsc[4,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,2,ic1,ic2,ic3,ic4]
               lhsc[4,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,3,ic1,ic2,ic3,ic4]
               lhsc[4,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,4,ic1,ic2,ic3,ic4]
               lhsc[4,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,5,ic1,ic2,ic3,ic4]
               rhs[4,ir1,ir2,ir3,ir4] -= coeff*rhs[2,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[5,2,i]
               lhsb[5,3,i] -= coeff*lhsb[2,3,i]
               lhsb[5,4,i] -= coeff*lhsb[2,4,i]
               lhsb[5,5,i] -= coeff*lhsb[2,5,i]
               lhsc[5,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,1,ic1,ic2,ic3,ic4]
               lhsc[5,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,2,ic1,ic2,ic3,ic4]
               lhsc[5,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,3,ic1,ic2,ic3,ic4]
               lhsc[5,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,4,ic1,ic2,ic3,ic4]
               lhsc[5,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,5,ic1,ic2,ic3,ic4]
               rhs[5,ir1,ir2,ir3,ir4] -= coeff*rhs[2,ir1,ir2,ir3,ir4]
         
               pivot = 1.00e0/lhsb[3,3,i]
               lhsb[3,4,i] = lhsb[3,4,i]*pivot
               lhsb[3,5,i] = lhsb[3,5,i]*pivot
               lhsc[3,1,ic1,ic2,ic3,ic4] *= pivot
               lhsc[3,2,ic1,ic2,ic3,ic4] *= pivot
               lhsc[3,3,ic1,ic2,ic3,ic4] *= pivot
               lhsc[3,4,ic1,ic2,ic3,ic4] *= pivot
               lhsc[3,5,ic1,ic2,ic3,ic4] *= pivot
               rhs[3,ir1,ir2,ir3,ir4] *= pivot
         
               coeff = lhsb[1,3,i]
               lhsb[1,4,i] -= coeff*lhsb[3,4,i]
               lhsb[1,5,i] -= coeff*lhsb[3,5,i]
               lhsc[1,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,1,ic1,ic2,ic3,ic4]
               lhsc[1,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,2,ic1,ic2,ic3,ic4]
               lhsc[1,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,3,ic1,ic2,ic3,ic4]
               lhsc[1,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,4,ic1,ic2,ic3,ic4]
               lhsc[1,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,5,ic1,ic2,ic3,ic4]
               rhs[1,ir1,ir2,ir3,ir4] -= coeff*rhs[3,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[2,3,i]
               lhsb[2,4,i] -= coeff*lhsb[3,4,i]
               lhsb[2,5,i] -= coeff*lhsb[3,5,i]
               lhsc[2,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,1,ic1,ic2,ic3,ic4]
               lhsc[2,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,2,ic1,ic2,ic3,ic4]
               lhsc[2,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,3,ic1,ic2,ic3,ic4]
               lhsc[2,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,4,ic1,ic2,ic3,ic4]
               lhsc[2,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,5,ic1,ic2,ic3,ic4]
               rhs[2,ir1,ir2,ir3,ir4] -= coeff*rhs[3,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[4,3,i]
               lhsb[4,4,i] -= coeff*lhsb[3,4,i]
               lhsb[4,5,i] -= coeff*lhsb[3,5,i]
               lhsc[4,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,1,ic1,ic2,ic3,ic4]
               lhsc[4,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,2,ic1,ic2,ic3,ic4]
               lhsc[4,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,3,ic1,ic2,ic3,ic4]
               lhsc[4,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,4,ic1,ic2,ic3,ic4]
               lhsc[4,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,5,ic1,ic2,ic3,ic4]
               rhs[4,ir1,ir2,ir3,ir4] -= coeff*rhs[3,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[5,3,i]
               lhsb[5,4,i] -= coeff*lhsb[3,4,i]
               lhsb[5,5,i] -= coeff*lhsb[3,5,i]
               lhsc[5,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,1,ic1,ic2,ic3,ic4]
               lhsc[5,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,2,ic1,ic2,ic3,ic4]
               lhsc[5,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,3,ic1,ic2,ic3,ic4]
               lhsc[5,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,4,ic1,ic2,ic3,ic4]
               lhsc[5,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,5,ic1,ic2,ic3,ic4]
               rhs[5,ir1,ir2,ir3,ir4] -= coeff*rhs[3,ir1,ir2,ir3,ir4]
         
               pivot = 1.00e0/lhsb[4,4,i]
               lhsb[4,5,i] = lhsb[4,5,i]*pivot
               lhsc[4,1,ic1,ic2,ic3,ic4] *= pivot
               lhsc[4,2,ic1,ic2,ic3,ic4] *= pivot
               lhsc[4,3,ic1,ic2,ic3,ic4] *= pivot
               lhsc[4,4,ic1,ic2,ic3,ic4] *= pivot
               lhsc[4,5,ic1,ic2,ic3,ic4] *= pivot
               rhs[4,ir1,ir2,ir3,ir4] *= pivot
         
               coeff = lhsb[1,4,i]
               lhsb[1,5,i] -= coeff*lhsb[4,5,i]
               lhsc[1,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,1,ic1,ic2,ic3,ic4]
               lhsc[1,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,2,ic1,ic2,ic3,ic4]
               lhsc[1,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,3,ic1,ic2,ic3,ic4]
               lhsc[1,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,4,ic1,ic2,ic3,ic4]
               lhsc[1,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,5,ic1,ic2,ic3,ic4]
               rhs[1,ir1,ir2,ir3,ir4] -= coeff*rhs[4,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[2,4,i]
               lhsb[2,5,i] -= coeff*lhsb[4,5,i]
               lhsc[2,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,1,ic1,ic2,ic3,ic4]
               lhsc[2,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,2,ic1,ic2,ic3,ic4]
               lhsc[2,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,3,ic1,ic2,ic3,ic4]
               lhsc[2,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,4,ic1,ic2,ic3,ic4]
               lhsc[2,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,5,ic1,ic2,ic3,ic4]
               rhs[2,ir1,ir2,ir3,ir4] -= coeff*rhs[4,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[3,4,i]
               lhsb[3,5,i] -= coeff*lhsb[4,5,i]
               lhsc[3,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,1,ic1,ic2,ic3,ic4]
               lhsc[3,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,2,ic1,ic2,ic3,ic4]
               lhsc[3,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,3,ic1,ic2,ic3,ic4]
               lhsc[3,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,4,ic1,ic2,ic3,ic4]
               lhsc[3,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,5,ic1,ic2,ic3,ic4]
               rhs[3,ir1,ir2,ir3,ir4] -= coeff*rhs[4,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[5,4,i]
               lhsb[5,5,i] -= coeff*lhsb[4,5,i]
               lhsc[5,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,1,ic1,ic2,ic3,ic4]
               lhsc[5,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,2,ic1,ic2,ic3,ic4]
               lhsc[5,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,3,ic1,ic2,ic3,ic4]
               lhsc[5,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,4,ic1,ic2,ic3,ic4]
               lhsc[5,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,5,ic1,ic2,ic3,ic4]
               rhs[5,ir1,ir2,ir3,ir4] -= coeff*rhs[4,ir1,ir2,ir3,ir4]
         
         
               pivot = 1.00e0/lhsb[5,5,i]
               lhsc[5,1,ic1,ic2,ic3,ic4] *= pivot
               lhsc[5,2,ic1,ic2,ic3,ic4] *= pivot
               lhsc[5,3,ic1,ic2,ic3,ic4] *= pivot
               lhsc[5,4,ic1,ic2,ic3,ic4] *= pivot
               lhsc[5,5,ic1,ic2,ic3,ic4] *= pivot
               rhs[5,ir1,ir2,ir3,ir4]   *= pivot
         
               coeff = lhsb[1,5,i]
               lhsc[1,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,1,ic1,ic2,ic3,ic4]
               lhsc[1,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,2,ic1,ic2,ic3,ic4]
               lhsc[1,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,3,ic1,ic2,ic3,ic4]
               lhsc[1,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,4,ic1,ic2,ic3,ic4]
               lhsc[1,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,5,ic1,ic2,ic3,ic4]
               rhs[1,ir1,ir2,ir3,ir4] -= coeff*rhs[5,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[2,5,i]
               lhsc[2,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,1,ic1,ic2,ic3,ic4]
               lhsc[2,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,2,ic1,ic2,ic3,ic4]
               lhsc[2,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,3,ic1,ic2,ic3,ic4]
               lhsc[2,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,4,ic1,ic2,ic3,ic4]
               lhsc[2,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,5,ic1,ic2,ic3,ic4]
               rhs[2,ir1,ir2,ir3,ir4] -= coeff*rhs[5,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[3,5,i]
               lhsc[3,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,1,ic1,ic2,ic3,ic4]
               lhsc[3,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,2,ic1,ic2,ic3,ic4]
               lhsc[3,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,3,ic1,ic2,ic3,ic4]
               lhsc[3,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,4,ic1,ic2,ic3,ic4]
               lhsc[3,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,5,ic1,ic2,ic3,ic4]
               rhs[3,ir1,ir2,ir3,ir4] -= coeff*rhs[5,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[4,5,i]
               lhsc[4,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,1,ic1,ic2,ic3,ic4]
               lhsc[4,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,2,ic1,ic2,ic3,ic4]
               lhsc[4,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,3,ic1,ic2,ic3,ic4]
               lhsc[4,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,4,ic1,ic2,ic3,ic4]
               lhsc[4,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,5,ic1,ic2,ic3,ic4]
               rhs[4,ir1,ir2,ir3,ir4] -= coeff*rhs[5,ir1,ir2,ir3,ir4] =#
         

            end

#---------------------------------------------------------------------
#     begin inner most do loop
#     do all the elements of the cell unless last 
#---------------------------------------------------------------------
            for i = istart+FIRST:isize-LAST

#---------------------------------------------------------------------
#     rhs[i] = rhs[i] - A*rhs[i-1]
#---------------------------------------------------------------------
#               matvec_sub(view(lhsa, 1:5, 1:5, i), view(rhs, 1:5, i-1, j, k, c), view(rhs, 1:5, i, j, k, c))
               matvec_sub(lhsa, i, rhs, i-1, j, k, c, rhs, i, j, k, c)

              #= a1, a2, a3, a4 = i-1, j, k, c
               b1, b2, b3, b4 = i, j, k, c

               rhs[1,b1,b2,b3,b4] -= lhsa[1,1,i]*rhs[1,a1,a2,a3,a4] +
                                       lhsa[1,2,i]*rhs[2,a1,a2,a3,a4] +
                                       lhsa[1,3,i]*rhs[3,a1,a2,a3,a4] +
                                       lhsa[1,4,i]*rhs[4,a1,a2,a3,a4] +
                                       lhsa[1,5,i]*rhs[5,a1,a2,a3,a4]
               rhs[2,b1,b2,b3,b4] -= lhsa[2,1,i]*rhs[1,a1,a2,a3,a4] +
                                       lhsa[2,2,i]*rhs[2,a1,a2,a3,a4] +
                                       lhsa[2,3,i]*rhs[3,a1,a2,a3,a4] +
                                       lhsa[2,4,i]*rhs[4,a1,a2,a3,a4] +
                                       lhsa[2,5,i]*rhs[5,a1,a2,a3,a4]
               rhs[3,b1,b2,b3,b4] -= lhsa[3,1,i]*rhs[1,a1,a2,a3,a4] +
                                       lhsa[3,2,i]*rhs[2,a1,a2,a3,a4] +
                                       lhsa[3,3,i]*rhs[3,a1,a2,a3,a4] +
                                       lhsa[3,4,i]*rhs[4,a1,a2,a3,a4] +
                                       lhsa[3,5,i]*rhs[5,a1,a2,a3,a4]
               rhs[4,b1,b2,b3,b4] -= lhsa[4,1,i]*rhs[1,a1,a2,a3,a4] +
                                       lhsa[4,2,i]*rhs[2,a1,a2,a3,a4] +
                                       lhsa[4,3,i]*rhs[3,a1,a2,a3,a4] +
                                       lhsa[4,4,i]*rhs[4,a1,a2,a3,a4] +
                                       lhsa[4,5,i]*rhs[5,a1,a2,a3,a4]
               rhs[5,b1,b2,b3,b4] -= lhsa[5,1,i]*rhs[1,a1,a2,a3,a4] +
                                       lhsa[5,2,i]*rhs[2,a1,a2,a3,a4] +
                                       lhsa[5,3,i]*rhs[3,a1,a2,a3,a4] +
                                       lhsa[5,4,i]*rhs[4,a1,a2,a3,a4] +
                                       lhsa[5,5,i]*rhs[5,a1,a2,a3,a4] =#


#---------------------------------------------------------------------
#      b[i] =  b[i] - C(i-1)*A(i)
#---------------------------------------------------------------------
               matmul_sub2(view(lhsa, 1:5, 1:5, i), view(lhsc, 1:5, 1:5, i-1, j, k, c), view(lhsb, 1:5, 1:5, i))
#               matmul_sub(lhsa, i, lhsc, i-1, j, k, c, lhsb, i)

             #=  ia = i 
               b1, b2, b3, b4 = i-1, j, k, c
               ic = i

               lhsb[1,1,ic] -= lhsa[1,1,ia]*lhsc[1,1,b1,b2,b3,b4] +
                                 lhsa[1,2,ia]*lhsc[2,1,b1,b2,b3,b4] +
                                 lhsa[1,3,ia]*lhsc[3,1,b1,b2,b3,b4] +
                                 lhsa[1,4,ia]*lhsc[4,1,b1,b2,b3,b4] +
                                 lhsa[1,5,ia]*lhsc[5,1,b1,b2,b3,b4]
               lhsb[2,1,ic] -= lhsa[2,1,ia]*lhsc[1,1,b1,b2,b3,b4] +
                                 lhsa[2,2,ia]*lhsc[2,1,b1,b2,b3,b4] +
                                 lhsa[2,3,ia]*lhsc[3,1,b1,b2,b3,b4] +
                                 lhsa[2,4,ia]*lhsc[4,1,b1,b2,b3,b4] +
                                 lhsa[2,5,ia]*lhsc[5,1,b1,b2,b3,b4]
               lhsb[3,1,ic] -= lhsa[3,1,ia]*lhsc[1,1,b1,b2,b3,b4] +
                                 lhsa[3,2,ia]*lhsc[2,1,b1,b2,b3,b4] +
                                 lhsa[3,3,ia]*lhsc[3,1,b1,b2,b3,b4] +
                                 lhsa[3,4,ia]*lhsc[4,1,b1,b2,b3,b4] +
                                 lhsa[3,5,ia]*lhsc[5,1,b1,b2,b3,b4]
               lhsb[4,1,ic] -= lhsa[4,1,ia]*lhsc[1,1,b1,b2,b3,b4] +
                                 lhsa[4,2,ia]*lhsc[2,1,b1,b2,b3,b4] +
                                 lhsa[4,3,ia]*lhsc[3,1,b1,b2,b3,b4] +
                                 lhsa[4,4,ia]*lhsc[4,1,b1,b2,b3,b4] +
                                 lhsa[4,5,ia]*lhsc[5,1,b1,b2,b3,b4]
               lhsb[5,1,ic] -= lhsa[5,1,ia]*lhsc[1,1,b1,b2,b3,b4] +
                                 lhsa[5,2,ia]*lhsc[2,1,b1,b2,b3,b4] +
                                 lhsa[5,3,ia]*lhsc[3,1,b1,b2,b3,b4] +
                                 lhsa[5,4,ia]*lhsc[4,1,b1,b2,b3,b4] +
                                 lhsa[5,5,ia]*lhsc[5,1,b1,b2,b3,b4]
               lhsb[1,2,ic] -= lhsa[1,1,ia]*lhsc[1,2,b1,b2,b3,b4] +
                                 lhsa[1,2,ia]*lhsc[2,2,b1,b2,b3,b4] +
                                 lhsa[1,3,ia]*lhsc[3,2,b1,b2,b3,b4] +
                                 lhsa[1,4,ia]*lhsc[4,2,b1,b2,b3,b4] +
                                 lhsa[1,5,ia]*lhsc[5,2,b1,b2,b3,b4]
               lhsb[2,2,ic] -= lhsa[2,1,ia]*lhsc[1,2,b1,b2,b3,b4] +
                                 lhsa[2,2,ia]*lhsc[2,2,b1,b2,b3,b4] +
                                 lhsa[2,3,ia]*lhsc[3,2,b1,b2,b3,b4] +
                                 lhsa[2,4,ia]*lhsc[4,2,b1,b2,b3,b4] +
                                 lhsa[2,5,ia]*lhsc[5,2,b1,b2,b3,b4]
               lhsb[3,2,ic] -= lhsa[3,1,ia]*lhsc[1,2,b1,b2,b3,b4] +
                                 lhsa[3,2,ia]*lhsc[2,2,b1,b2,b3,b4] +
                                 lhsa[3,3,ia]*lhsc[3,2,b1,b2,b3,b4] +
                                 lhsa[3,4,ia]*lhsc[4,2,b1,b2,b3,b4] +
                                 lhsa[3,5,ia]*lhsc[5,2,b1,b2,b3,b4]
               lhsb[4,2,ic] -= lhsa[4,1,ia]*lhsc[1,2,b1,b2,b3,b4] +
                                 lhsa[4,2,ia]*lhsc[2,2,b1,b2,b3,b4] +
                                 lhsa[4,3,ia]*lhsc[3,2,b1,b2,b3,b4] +
                                 lhsa[4,4,ia]*lhsc[4,2,b1,b2,b3,b4] +
                                 lhsa[4,5,ia]*lhsc[5,2,b1,b2,b3,b4]
               lhsb[5,2,ic] -= lhsa[5,1,ia]*lhsc[1,2,b1,b2,b3,b4] +
                                 lhsa[5,2,ia]*lhsc[2,2,b1,b2,b3,b4] +
                                 lhsa[5,3,ia]*lhsc[3,2,b1,b2,b3,b4] +
                                 lhsa[5,4,ia]*lhsc[4,2,b1,b2,b3,b4] +
                                 lhsa[5,5,ia]*lhsc[5,2,b1,b2,b3,b4]
               lhsb[1,3,ic] -= lhsa[1,1,ia]*lhsc[1,3,b1,b2,b3,b4] +
                                 lhsa[1,2,ia]*lhsc[2,3,b1,b2,b3,b4] +
                                 lhsa[1,3,ia]*lhsc[3,3,b1,b2,b3,b4] +
                                 lhsa[1,4,ia]*lhsc[4,3,b1,b2,b3,b4] +
                                 lhsa[1,5,ia]*lhsc[5,3,b1,b2,b3,b4]
               lhsb[2,3,ic] -= lhsa[2,1,ia]*lhsc[1,3,b1,b2,b3,b4] +
                                 lhsa[2,2,ia]*lhsc[2,3,b1,b2,b3,b4] +
                                 lhsa[2,3,ia]*lhsc[3,3,b1,b2,b3,b4] +
                                 lhsa[2,4,ia]*lhsc[4,3,b1,b2,b3,b4] +
                                 lhsa[2,5,ia]*lhsc[5,3,b1,b2,b3,b4]
               lhsb[3,3,ic] -= lhsa[3,1,ia]*lhsc[1,3,b1,b2,b3,b4] +
                                 lhsa[3,2,ia]*lhsc[2,3,b1,b2,b3,b4] +
                                 lhsa[3,3,ia]*lhsc[3,3,b1,b2,b3,b4] +
                                 lhsa[3,4,ia]*lhsc[4,3,b1,b2,b3,b4] +
                                 lhsa[3,5,ia]*lhsc[5,3,b1,b2,b3,b4]
               lhsb[4,3,ic] -= lhsa[4,1,ia]*lhsc[1,3,b1,b2,b3,b4] +
                                 lhsa[4,2,ia]*lhsc[2,3,b1,b2,b3,b4] +
                                 lhsa[4,3,ia]*lhsc[3,3,b1,b2,b3,b4] +
                                 lhsa[4,4,ia]*lhsc[4,3,b1,b2,b3,b4] +
                                 lhsa[4,5,ia]*lhsc[5,3,b1,b2,b3,b4]
               lhsb[5,3,ic] -= lhsa[5,1,ia]*lhsc[1,3,b1,b2,b3,b4] +
                                 lhsa[5,2,ia]*lhsc[2,3,b1,b2,b3,b4] +
                                 lhsa[5,3,ia]*lhsc[3,3,b1,b2,b3,b4] +
                                 lhsa[5,4,ia]*lhsc[4,3,b1,b2,b3,b4] +
                                 lhsa[5,5,ia]*lhsc[5,3,b1,b2,b3,b4]
               lhsb[1,4,ic] -= lhsa[1,1,ia]*lhsc[1,4,b1,b2,b3,b4] +
                                 lhsa[1,2,ia]*lhsc[2,4,b1,b2,b3,b4] +
                                 lhsa[1,3,ia]*lhsc[3,4,b1,b2,b3,b4] +
                                 lhsa[1,4,ia]*lhsc[4,4,b1,b2,b3,b4] +
                                 lhsa[1,5,ia]*lhsc[5,4,b1,b2,b3,b4]
               lhsb[2,4,ic] -= lhsa[2,1,ia]*lhsc[1,4,b1,b2,b3,b4] +
                                 lhsa[2,2,ia]*lhsc[2,4,b1,b2,b3,b4] +
                                 lhsa[2,3,ia]*lhsc[3,4,b1,b2,b3,b4] +
                                 lhsa[2,4,ia]*lhsc[4,4,b1,b2,b3,b4] +
                                 lhsa[2,5,ia]*lhsc[5,4,b1,b2,b3,b4]
               lhsb[3,4,ic] -= lhsa[3,1,ia]*lhsc[1,4,b1,b2,b3,b4] +
                                 lhsa[3,2,ia]*lhsc[2,4,b1,b2,b3,b4] +
                                 lhsa[3,3,ia]*lhsc[3,4,b1,b2,b3,b4] +
                                 lhsa[3,4,ia]*lhsc[4,4,b1,b2,b3,b4] +
                                 lhsa[3,5,ia]*lhsc[5,4,b1,b2,b3,b4]
               lhsb[4,4,ic] -= lhsa[4,1,ia]*lhsc[1,4,b1,b2,b3,b4] +
                                 lhsa[4,2,ia]*lhsc[2,4,b1,b2,b3,b4] +
                                 lhsa[4,3,ia]*lhsc[3,4,b1,b2,b3,b4] +
                                 lhsa[4,4,ia]*lhsc[4,4,b1,b2,b3,b4] +
                                 lhsa[4,5,ia]*lhsc[5,4,b1,b2,b3,b4]
               lhsb[5,4,ic] -= lhsa[5,1,ia]*lhsc[1,4,b1,b2,b3,b4] +
                                 lhsa[5,2,ia]*lhsc[2,4,b1,b2,b3,b4] +
                                 lhsa[5,3,ia]*lhsc[3,4,b1,b2,b3,b4] +
                                 lhsa[5,4,ia]*lhsc[4,4,b1,b2,b3,b4] +
                                 lhsa[5,5,ia]*lhsc[5,4,b1,b2,b3,b4]
               lhsb[1,5,ic] -= lhsa[1,1,ia]*lhsc[1,5,b1,b2,b3,b4] +
                                 lhsa[1,2,ia]*lhsc[2,5,b1,b2,b3,b4] +
                                 lhsa[1,3,ia]*lhsc[3,5,b1,b2,b3,b4] +
                                 lhsa[1,4,ia]*lhsc[4,5,b1,b2,b3,b4] +
                                 lhsa[1,5,ia]*lhsc[5,5,b1,b2,b3,b4]
               lhsb[2,5,ic] -= lhsa[2,1,ia]*lhsc[1,5,b1,b2,b3,b4] +
                                 lhsa[2,2,ia]*lhsc[2,5,b1,b2,b3,b4] +
                                 lhsa[2,3,ia]*lhsc[3,5,b1,b2,b3,b4] +
                                 lhsa[2,4,ia]*lhsc[4,5,b1,b2,b3,b4] +
                                 lhsa[2,5,ia]*lhsc[5,5,b1,b2,b3,b4]
               lhsb[3,5,ic] -= lhsa[3,1,ia]*lhsc[1,5,b1,b2,b3,b4] +
                                 lhsa[3,2,ia]*lhsc[2,5,b1,b2,b3,b4] +
                                 lhsa[3,3,ia]*lhsc[3,5,b1,b2,b3,b4] +
                                 lhsa[3,4,ia]*lhsc[4,5,b1,b2,b3,b4] +
                                 lhsa[3,5,ia]*lhsc[5,5,b1,b2,b3,b4]
               lhsb[4,5,ic] -= lhsa[4,1,ia]*lhsc[1,5,b1,b2,b3,b4] +
                                 lhsa[4,2,ia]*lhsc[2,5,b1,b2,b3,b4] +
                                 lhsa[4,3,ia]*lhsc[3,5,b1,b2,b3,b4] +
                                 lhsa[4,4,ia]*lhsc[4,5,b1,b2,b3,b4] +
                                 lhsa[4,5,ia]*lhsc[5,5,b1,b2,b3,b4]
               lhsb[5,5,ic] -= lhsa[5,1,ia]*lhsc[1,5,b1,b2,b3,b4] +
                                 lhsa[5,2,ia]*lhsc[2,5,b1,b2,b3,b4] +
                                 lhsa[5,3,ia]*lhsc[3,5,b1,b2,b3,b4] +
                                 lhsa[5,4,ia]*lhsc[4,5,b1,b2,b3,b4] +
                                 lhsa[5,5,ia]*lhsc[5,5,b1,b2,b3,b4] =#


#---------------------------------------------------------------------
#     multiply c[i,j,k] by b_inverse and copy back to c
#     multiply rhs[1,j,k] by b_inverse(1,j,k) and copy to rhs
#---------------------------------------------------------------------
#               binvcrhs(view(lhsb, 1:5, 1:5, i), view(lhsc, 1:5, 1:5, i, j, k, c), view(rhs, 1:5, i, j, k, c))
               binvcrhs(lhsb, i, lhsc, i, j, k, c, rhs, i, j, k, c)

              #= ic1, ic2, ic3, ic4 = i, j, k, c
               ir1, ir2, ir3, ir4 = i, j, k, c

               pivot = 1.00e0/lhsb[1,1,i]
               lhsb[1,2,i] *= pivot
               lhsb[1,3,i] *= pivot
               lhsb[1,4,i] *= pivot
               lhsb[1,5,i] *= pivot
               lhsc[1,1,ic1,ic2,ic3,ic4] *= pivot
               lhsc[1,2,ic1,ic2,ic3,ic4] *= pivot
               lhsc[1,3,ic1,ic2,ic3,ic4] *= pivot
               lhsc[1,4,ic1,ic2,ic3,ic4] *= pivot
               lhsc[1,5,ic1,ic2,ic3,ic4] *= pivot
               rhs[1,ir1,ir2,ir3,ir4]   *= pivot
         
               coeff = lhsb[2,1,i]
               lhsb[2,2,i] -= coeff*lhsb[1,2,i]
               lhsb[2,3,i] -= coeff*lhsb[1,3,i]
               lhsb[2,4,i] -= coeff*lhsb[1,4,i]
               lhsb[2,5,i] -= coeff*lhsb[1,5,i]
               lhsc[2,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,1,ic1,ic2,ic3,ic4]
               lhsc[2,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,2,ic1,ic2,ic3,ic4]
               lhsc[2,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,3,ic1,ic2,ic3,ic4]
               lhsc[2,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,4,ic1,ic2,ic3,ic4]
               lhsc[2,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,5,ic1,ic2,ic3,ic4]
               rhs[2,ir1,ir2,ir3,ir4] -= coeff*rhs[1,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[3,1,i]
               lhsb[3,2,i] -= coeff*lhsb[1,2,i]
               lhsb[3,3,i] -= coeff*lhsb[1,3,i]
               lhsb[3,4,i] -= coeff*lhsb[1,4,i]
               lhsb[3,5,i] -= coeff*lhsb[1,5,i]
               lhsc[3,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,1,ic1,ic2,ic3,ic4]
               lhsc[3,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,2,ic1,ic2,ic3,ic4]
               lhsc[3,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,3,ic1,ic2,ic3,ic4]
               lhsc[3,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,4,ic1,ic2,ic3,ic4]
               lhsc[3,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,5,ic1,ic2,ic3,ic4]
               rhs[3,ir1,ir2,ir3,ir4] -= coeff*rhs[1,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[4,1,i]
               lhsb[4,2,i] -= coeff*lhsb[1,2,i]
               lhsb[4,3,i] -= coeff*lhsb[1,3,i]
               lhsb[4,4,i] -= coeff*lhsb[1,4,i]
               lhsb[4,5,i] -= coeff*lhsb[1,5,i]
               lhsc[4,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,1,ic1,ic2,ic3,ic4]
               lhsc[4,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,2,ic1,ic2,ic3,ic4]
               lhsc[4,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,3,ic1,ic2,ic3,ic4]
               lhsc[4,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,4,ic1,ic2,ic3,ic4]
               lhsc[4,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,5,ic1,ic2,ic3,ic4]
               rhs[4,ir1,ir2,ir3,ir4] -= coeff*rhs[1,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[5,1,i]
               lhsb[5,2,i] -= coeff*lhsb[1,2,i]
               lhsb[5,3,i] -= coeff*lhsb[1,3,i]
               lhsb[5,4,i] -= coeff*lhsb[1,4,i]
               lhsb[5,5,i] -= coeff*lhsb[1,5,i]
               lhsc[5,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,1,ic1,ic2,ic3,ic4]
               lhsc[5,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,2,ic1,ic2,ic3,ic4]
               lhsc[5,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,3,ic1,ic2,ic3,ic4]
               lhsc[5,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,4,ic1,ic2,ic3,ic4]
               lhsc[5,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[1,5,ic1,ic2,ic3,ic4]
               rhs[5,ir1,ir2,ir3,ir4] -= coeff*rhs[1,ir1,ir2,ir3,ir4]
         
         
               pivot = 1.00e0/lhsb[2,2,i]
               lhsb[2,3,i] *= pivot
               lhsb[2,4,i] *= pivot
               lhsb[2,5,i] *= pivot
               lhsc[2,1,ic1,ic2,ic3,ic4] *= pivot
               lhsc[2,2,ic1,ic2,ic3,ic4] *= pivot
               lhsc[2,3,ic1,ic2,ic3,ic4] *= pivot
               lhsc[2,4,ic1,ic2,ic3,ic4] *= pivot
               lhsc[2,5,ic1,ic2,ic3,ic4] *= pivot
               rhs[2,ir1,ir2,ir3,ir4] *= pivot
         
               coeff = lhsb[1,2,i]
               lhsb[1,3,i] -= coeff*lhsb[2,3,i]
               lhsb[1,4,i] -= coeff*lhsb[2,4,i]
               lhsb[1,5,i] -= coeff*lhsb[2,5,i]
               lhsc[1,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,1,ic1,ic2,ic3,ic4]
               lhsc[1,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,2,ic1,ic2,ic3,ic4]
               lhsc[1,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,3,ic1,ic2,ic3,ic4]
               lhsc[1,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,4,ic1,ic2,ic3,ic4]
               lhsc[1,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,5,ic1,ic2,ic3,ic4]
               rhs[1,ir1,ir2,ir3,ir4] -= coeff*rhs[2,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[3,2,i]
               lhsb[3,3,i] -= coeff*lhsb[2,3,i]
               lhsb[3,4,i] -= coeff*lhsb[2,4,i]
               lhsb[3,5,i] -= coeff*lhsb[2,5,i]
               lhsc[3,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,1,ic1,ic2,ic3,ic4]
               lhsc[3,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,2,ic1,ic2,ic3,ic4]
               lhsc[3,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,3,ic1,ic2,ic3,ic4]
               lhsc[3,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,4,ic1,ic2,ic3,ic4]
               lhsc[3,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,5,ic1,ic2,ic3,ic4]
               rhs[3,ir1,ir2,ir3,ir4] -= coeff*rhs[2,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[4,2,i]
               lhsb[4,3,i] -= coeff*lhsb[2,3,i]
               lhsb[4,4,i] -= coeff*lhsb[2,4,i]
               lhsb[4,5,i] -= coeff*lhsb[2,5,i]
               lhsc[4,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,1,ic1,ic2,ic3,ic4]
               lhsc[4,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,2,ic1,ic2,ic3,ic4]
               lhsc[4,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,3,ic1,ic2,ic3,ic4]
               lhsc[4,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,4,ic1,ic2,ic3,ic4]
               lhsc[4,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,5,ic1,ic2,ic3,ic4]
               rhs[4,ir1,ir2,ir3,ir4] -= coeff*rhs[2,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[5,2,i]
               lhsb[5,3,i] -= coeff*lhsb[2,3,i]
               lhsb[5,4,i] -= coeff*lhsb[2,4,i]
               lhsb[5,5,i] -= coeff*lhsb[2,5,i]
               lhsc[5,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,1,ic1,ic2,ic3,ic4]
               lhsc[5,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,2,ic1,ic2,ic3,ic4]
               lhsc[5,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,3,ic1,ic2,ic3,ic4]
               lhsc[5,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,4,ic1,ic2,ic3,ic4]
               lhsc[5,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[2,5,ic1,ic2,ic3,ic4]
               rhs[5,ir1,ir2,ir3,ir4] -= coeff*rhs[2,ir1,ir2,ir3,ir4]
         
               pivot = 1.00e0/lhsb[3,3,i]
               lhsb[3,4,i] = lhsb[3,4,i]*pivot
               lhsb[3,5,i] = lhsb[3,5,i]*pivot
               lhsc[3,1,ic1,ic2,ic3,ic4] *= pivot
               lhsc[3,2,ic1,ic2,ic3,ic4] *= pivot
               lhsc[3,3,ic1,ic2,ic3,ic4] *= pivot
               lhsc[3,4,ic1,ic2,ic3,ic4] *= pivot
               lhsc[3,5,ic1,ic2,ic3,ic4] *= pivot
               rhs[3,ir1,ir2,ir3,ir4] *= pivot
         
               coeff = lhsb[1,3,i]
               lhsb[1,4,i] -= coeff*lhsb[3,4,i]
               lhsb[1,5,i] -= coeff*lhsb[3,5,i]
               lhsc[1,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,1,ic1,ic2,ic3,ic4]
               lhsc[1,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,2,ic1,ic2,ic3,ic4]
               lhsc[1,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,3,ic1,ic2,ic3,ic4]
               lhsc[1,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,4,ic1,ic2,ic3,ic4]
               lhsc[1,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,5,ic1,ic2,ic3,ic4]
               rhs[1,ir1,ir2,ir3,ir4] -= coeff*rhs[3,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[2,3,i]
               lhsb[2,4,i] -= coeff*lhsb[3,4,i]
               lhsb[2,5,i] -= coeff*lhsb[3,5,i]
               lhsc[2,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,1,ic1,ic2,ic3,ic4]
               lhsc[2,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,2,ic1,ic2,ic3,ic4]
               lhsc[2,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,3,ic1,ic2,ic3,ic4]
               lhsc[2,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,4,ic1,ic2,ic3,ic4]
               lhsc[2,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,5,ic1,ic2,ic3,ic4]
               rhs[2,ir1,ir2,ir3,ir4] -= coeff*rhs[3,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[4,3,i]
               lhsb[4,4,i] -= coeff*lhsb[3,4,i]
               lhsb[4,5,i] -= coeff*lhsb[3,5,i]
               lhsc[4,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,1,ic1,ic2,ic3,ic4]
               lhsc[4,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,2,ic1,ic2,ic3,ic4]
               lhsc[4,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,3,ic1,ic2,ic3,ic4]
               lhsc[4,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,4,ic1,ic2,ic3,ic4]
               lhsc[4,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,5,ic1,ic2,ic3,ic4]
               rhs[4,ir1,ir2,ir3,ir4] -= coeff*rhs[3,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[5,3,i]
               lhsb[5,4,i] -= coeff*lhsb[3,4,i]
               lhsb[5,5,i] -= coeff*lhsb[3,5,i]
               lhsc[5,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,1,ic1,ic2,ic3,ic4]
               lhsc[5,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,2,ic1,ic2,ic3,ic4]
               lhsc[5,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,3,ic1,ic2,ic3,ic4]
               lhsc[5,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,4,ic1,ic2,ic3,ic4]
               lhsc[5,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[3,5,ic1,ic2,ic3,ic4]
               rhs[5,ir1,ir2,ir3,ir4] -= coeff*rhs[3,ir1,ir2,ir3,ir4]
         
               pivot = 1.00e0/lhsb[4,4,i]
               lhsb[4,5,i] = lhsb[4,5,i]*pivot
               lhsc[4,1,ic1,ic2,ic3,ic4] *= pivot
               lhsc[4,2,ic1,ic2,ic3,ic4] *= pivot
               lhsc[4,3,ic1,ic2,ic3,ic4] *= pivot
               lhsc[4,4,ic1,ic2,ic3,ic4] *= pivot
               lhsc[4,5,ic1,ic2,ic3,ic4] *= pivot
               rhs[4,ir1,ir2,ir3,ir4] *= pivot
         
               coeff = lhsb[1,4,i]
               lhsb[1,5,i] -= coeff*lhsb[4,5,i]
               lhsc[1,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,1,ic1,ic2,ic3,ic4]
               lhsc[1,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,2,ic1,ic2,ic3,ic4]
               lhsc[1,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,3,ic1,ic2,ic3,ic4]
               lhsc[1,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,4,ic1,ic2,ic3,ic4]
               lhsc[1,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,5,ic1,ic2,ic3,ic4]
               rhs[1,ir1,ir2,ir3,ir4] -= coeff*rhs[4,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[2,4,i]
               lhsb[2,5,i] -= coeff*lhsb[4,5,i]
               lhsc[2,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,1,ic1,ic2,ic3,ic4]
               lhsc[2,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,2,ic1,ic2,ic3,ic4]
               lhsc[2,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,3,ic1,ic2,ic3,ic4]
               lhsc[2,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,4,ic1,ic2,ic3,ic4]
               lhsc[2,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,5,ic1,ic2,ic3,ic4]
               rhs[2,ir1,ir2,ir3,ir4] -= coeff*rhs[4,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[3,4,i]
               lhsb[3,5,i] -= coeff*lhsb[4,5,i]
               lhsc[3,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,1,ic1,ic2,ic3,ic4]
               lhsc[3,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,2,ic1,ic2,ic3,ic4]
               lhsc[3,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,3,ic1,ic2,ic3,ic4]
               lhsc[3,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,4,ic1,ic2,ic3,ic4]
               lhsc[3,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,5,ic1,ic2,ic3,ic4]
               rhs[3,ir1,ir2,ir3,ir4] -= coeff*rhs[4,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[5,4,i]
               lhsb[5,5,i] -= coeff*lhsb[4,5,i]
               lhsc[5,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,1,ic1,ic2,ic3,ic4]
               lhsc[5,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,2,ic1,ic2,ic3,ic4]
               lhsc[5,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,3,ic1,ic2,ic3,ic4]
               lhsc[5,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,4,ic1,ic2,ic3,ic4]
               lhsc[5,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[4,5,ic1,ic2,ic3,ic4]
               rhs[5,ir1,ir2,ir3,ir4] -= coeff*rhs[4,ir1,ir2,ir3,ir4]
         
         
               pivot = 1.00e0/lhsb[5,5,i]
               lhsc[5,1,ic1,ic2,ic3,ic4] *= pivot
               lhsc[5,2,ic1,ic2,ic3,ic4] *= pivot
               lhsc[5,3,ic1,ic2,ic3,ic4] *= pivot
               lhsc[5,4,ic1,ic2,ic3,ic4] *= pivot
               lhsc[5,5,ic1,ic2,ic3,ic4] *= pivot
               rhs[5,ir1,ir2,ir3,ir4]   *= pivot
         
               coeff = lhsb[1,5,i]
               lhsc[1,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,1,ic1,ic2,ic3,ic4]
               lhsc[1,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,2,ic1,ic2,ic3,ic4]
               lhsc[1,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,3,ic1,ic2,ic3,ic4]
               lhsc[1,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,4,ic1,ic2,ic3,ic4]
               lhsc[1,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,5,ic1,ic2,ic3,ic4]
               rhs[1,ir1,ir2,ir3,ir4] -= coeff*rhs[5,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[2,5,i]
               lhsc[2,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,1,ic1,ic2,ic3,ic4]
               lhsc[2,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,2,ic1,ic2,ic3,ic4]
               lhsc[2,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,3,ic1,ic2,ic3,ic4]
               lhsc[2,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,4,ic1,ic2,ic3,ic4]
               lhsc[2,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,5,ic1,ic2,ic3,ic4]
               rhs[2,ir1,ir2,ir3,ir4] -= coeff*rhs[5,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[3,5,i]
               lhsc[3,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,1,ic1,ic2,ic3,ic4]
               lhsc[3,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,2,ic1,ic2,ic3,ic4]
               lhsc[3,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,3,ic1,ic2,ic3,ic4]
               lhsc[3,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,4,ic1,ic2,ic3,ic4]
               lhsc[3,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,5,ic1,ic2,ic3,ic4]
               rhs[3,ir1,ir2,ir3,ir4] -= coeff*rhs[5,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[4,5,i]
               lhsc[4,1,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,1,ic1,ic2,ic3,ic4]
               lhsc[4,2,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,2,ic1,ic2,ic3,ic4]
               lhsc[4,3,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,3,ic1,ic2,ic3,ic4]
               lhsc[4,4,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,4,ic1,ic2,ic3,ic4]
               lhsc[4,5,ic1,ic2,ic3,ic4] -= coeff*lhsc[5,5,ic1,ic2,ic3,ic4]
               rhs[4,ir1,ir2,ir3,ir4] -= coeff*rhs[5,ir1,ir2,ir3,ir4] =#
            end

#---------------------------------------------------------------------
#     Now finish up special cases for last cell
#---------------------------------------------------------------------
            if LAST == 1

#---------------------------------------------------------------------
#     rhs[isize] = rhs[isize] - A*rhs[isize-1]
#---------------------------------------------------------------------
#               matvec_sub(view(lhsa, 1:5, 1:5, isize), view(rhs, 1:5, isize-1, j, k, c), view(rhs, 1:5, isize, j, k, c))
               matvec_sub(lhsa, isize, rhs, isize-1, j, k, c, rhs, isize, j, k, c)

              #= i = isize
               a1, a2, a3, a4 = isize-1, j, k, c
               b1, b2, b3, b4 = isize, j, k, c

               rhs[1,b1,b2,b3,b4] -= lhsa[1,1,i]*rhs[1,a1,a2,a3,a4] +
                                       lhsa[1,2,i]*rhs[2,a1,a2,a3,a4] +
                                       lhsa[1,3,i]*rhs[3,a1,a2,a3,a4] +
                                       lhsa[1,4,i]*rhs[4,a1,a2,a3,a4] +
                                       lhsa[1,5,i]*rhs[5,a1,a2,a3,a4]
               rhs[2,b1,b2,b3,b4] -= lhsa[2,1,i]*rhs[1,a1,a2,a3,a4] +
                                       lhsa[2,2,i]*rhs[2,a1,a2,a3,a4] +
                                       lhsa[2,3,i]*rhs[3,a1,a2,a3,a4] +
                                       lhsa[2,4,i]*rhs[4,a1,a2,a3,a4] +
                                       lhsa[2,5,i]*rhs[5,a1,a2,a3,a4]
               rhs[3,b1,b2,b3,b4] -= lhsa[3,1,i]*rhs[1,a1,a2,a3,a4] +
                                       lhsa[3,2,i]*rhs[2,a1,a2,a3,a4] +
                                       lhsa[3,3,i]*rhs[3,a1,a2,a3,a4] +
                                       lhsa[3,4,i]*rhs[4,a1,a2,a3,a4] +
                                       lhsa[3,5,i]*rhs[5,a1,a2,a3,a4]
               rhs[4,b1,b2,b3,b4] -= lhsa[4,1,i]*rhs[1,a1,a2,a3,a4] +
                                       lhsa[4,2,i]*rhs[2,a1,a2,a3,a4] +
                                       lhsa[4,3,i]*rhs[3,a1,a2,a3,a4] +
                                       lhsa[4,4,i]*rhs[4,a1,a2,a3,a4] +
                                       lhsa[4,5,i]*rhs[5,a1,a2,a3,a4]
               rhs[5,b1,b2,b3,b4] -= lhsa[5,1,i]*rhs[1,a1,a2,a3,a4] +
                                       lhsa[5,2,i]*rhs[2,a1,a2,a3,a4] +
                                       lhsa[5,3,i]*rhs[3,a1,a2,a3,a4] +
                                       lhsa[5,4,i]*rhs[4,a1,a2,a3,a4] +
                                       lhsa[5,5,i]*rhs[5,a1,a2,a3,a4] =#


#---------------------------------------------------------------------
#      b[isize] =  b[isize] - C(isize-1)*A(isize)
#---------------------------------------------------------------------
               matmul_sub2(view(lhsa, 1:5, 1:5, isize), view(lhsc, 1:5, 1:5, isize-1, j, k, c), view(lhsb, 1:5, 1:5, isize))
#               matmul_sub(lhsa, isize, lhsc, isize-1, j, k, c, lhsb, isize)

               #=ia = i 
               b1, b2, b3, b4 = i-1, j, k, c
               ic = i

               lhsb[1,1,ic] -= lhsa[1,1,ia]*lhsc[1,1,b1,b2,b3,b4] +
                                 lhsa[1,2,ia]*lhsc[2,1,b1,b2,b3,b4] +
                                 lhsa[1,3,ia]*lhsc[3,1,b1,b2,b3,b4] +
                                 lhsa[1,4,ia]*lhsc[4,1,b1,b2,b3,b4] +
                                 lhsa[1,5,ia]*lhsc[5,1,b1,b2,b3,b4]
               lhsb[2,1,ic] -= lhsa[2,1,ia]*lhsc[1,1,b1,b2,b3,b4] +
                                 lhsa[2,2,ia]*lhsc[2,1,b1,b2,b3,b4] +
                                 lhsa[2,3,ia]*lhsc[3,1,b1,b2,b3,b4] +
                                 lhsa[2,4,ia]*lhsc[4,1,b1,b2,b3,b4] +
                                 lhsa[2,5,ia]*lhsc[5,1,b1,b2,b3,b4]
               lhsb[3,1,ic] -= lhsa[3,1,ia]*lhsc[1,1,b1,b2,b3,b4] +
                                 lhsa[3,2,ia]*lhsc[2,1,b1,b2,b3,b4] +
                                 lhsa[3,3,ia]*lhsc[3,1,b1,b2,b3,b4] +
                                 lhsa[3,4,ia]*lhsc[4,1,b1,b2,b3,b4] +
                                 lhsa[3,5,ia]*lhsc[5,1,b1,b2,b3,b4]
               lhsb[4,1,ic] -= lhsa[4,1,ia]*lhsc[1,1,b1,b2,b3,b4] +
                                 lhsa[4,2,ia]*lhsc[2,1,b1,b2,b3,b4] +
                                 lhsa[4,3,ia]*lhsc[3,1,b1,b2,b3,b4] +
                                 lhsa[4,4,ia]*lhsc[4,1,b1,b2,b3,b4] +
                                 lhsa[4,5,ia]*lhsc[5,1,b1,b2,b3,b4]
               lhsb[5,1,ic] -= lhsa[5,1,ia]*lhsc[1,1,b1,b2,b3,b4] +
                                 lhsa[5,2,ia]*lhsc[2,1,b1,b2,b3,b4] +
                                 lhsa[5,3,ia]*lhsc[3,1,b1,b2,b3,b4] +
                                 lhsa[5,4,ia]*lhsc[4,1,b1,b2,b3,b4] +
                                 lhsa[5,5,ia]*lhsc[5,1,b1,b2,b3,b4]
               lhsb[1,2,ic] -= lhsa[1,1,ia]*lhsc[1,2,b1,b2,b3,b4] +
                                 lhsa[1,2,ia]*lhsc[2,2,b1,b2,b3,b4] +
                                 lhsa[1,3,ia]*lhsc[3,2,b1,b2,b3,b4] +
                                 lhsa[1,4,ia]*lhsc[4,2,b1,b2,b3,b4] +
                                 lhsa[1,5,ia]*lhsc[5,2,b1,b2,b3,b4]
               lhsb[2,2,ic] -= lhsa[2,1,ia]*lhsc[1,2,b1,b2,b3,b4] +
                                 lhsa[2,2,ia]*lhsc[2,2,b1,b2,b3,b4] +
                                 lhsa[2,3,ia]*lhsc[3,2,b1,b2,b3,b4] +
                                 lhsa[2,4,ia]*lhsc[4,2,b1,b2,b3,b4] +
                                 lhsa[2,5,ia]*lhsc[5,2,b1,b2,b3,b4]
               lhsb[3,2,ic] -= lhsa[3,1,ia]*lhsc[1,2,b1,b2,b3,b4] +
                                 lhsa[3,2,ia]*lhsc[2,2,b1,b2,b3,b4] +
                                 lhsa[3,3,ia]*lhsc[3,2,b1,b2,b3,b4] +
                                 lhsa[3,4,ia]*lhsc[4,2,b1,b2,b3,b4] +
                                 lhsa[3,5,ia]*lhsc[5,2,b1,b2,b3,b4]
               lhsb[4,2,ic] -= lhsa[4,1,ia]*lhsc[1,2,b1,b2,b3,b4] +
                                 lhsa[4,2,ia]*lhsc[2,2,b1,b2,b3,b4] +
                                 lhsa[4,3,ia]*lhsc[3,2,b1,b2,b3,b4] +
                                 lhsa[4,4,ia]*lhsc[4,2,b1,b2,b3,b4] +
                                 lhsa[4,5,ia]*lhsc[5,2,b1,b2,b3,b4]
               lhsb[5,2,ic] -= lhsa[5,1,ia]*lhsc[1,2,b1,b2,b3,b4] +
                                 lhsa[5,2,ia]*lhsc[2,2,b1,b2,b3,b4] +
                                 lhsa[5,3,ia]*lhsc[3,2,b1,b2,b3,b4] +
                                 lhsa[5,4,ia]*lhsc[4,2,b1,b2,b3,b4] +
                                 lhsa[5,5,ia]*lhsc[5,2,b1,b2,b3,b4]
               lhsb[1,3,ic] -= lhsa[1,1,ia]*lhsc[1,3,b1,b2,b3,b4] +
                                 lhsa[1,2,ia]*lhsc[2,3,b1,b2,b3,b4] +
                                 lhsa[1,3,ia]*lhsc[3,3,b1,b2,b3,b4] +
                                 lhsa[1,4,ia]*lhsc[4,3,b1,b2,b3,b4] +
                                 lhsa[1,5,ia]*lhsc[5,3,b1,b2,b3,b4]
               lhsb[2,3,ic] -= lhsa[2,1,ia]*lhsc[1,3,b1,b2,b3,b4] +
                                 lhsa[2,2,ia]*lhsc[2,3,b1,b2,b3,b4] +
                                 lhsa[2,3,ia]*lhsc[3,3,b1,b2,b3,b4] +
                                 lhsa[2,4,ia]*lhsc[4,3,b1,b2,b3,b4] +
                                 lhsa[2,5,ia]*lhsc[5,3,b1,b2,b3,b4]
               lhsb[3,3,ic] -= lhsa[3,1,ia]*lhsc[1,3,b1,b2,b3,b4] +
                                 lhsa[3,2,ia]*lhsc[2,3,b1,b2,b3,b4] +
                                 lhsa[3,3,ia]*lhsc[3,3,b1,b2,b3,b4] +
                                 lhsa[3,4,ia]*lhsc[4,3,b1,b2,b3,b4] +
                                 lhsa[3,5,ia]*lhsc[5,3,b1,b2,b3,b4]
               lhsb[4,3,ic] -= lhsa[4,1,ia]*lhsc[1,3,b1,b2,b3,b4] +
                                 lhsa[4,2,ia]*lhsc[2,3,b1,b2,b3,b4] +
                                 lhsa[4,3,ia]*lhsc[3,3,b1,b2,b3,b4] +
                                 lhsa[4,4,ia]*lhsc[4,3,b1,b2,b3,b4] +
                                 lhsa[4,5,ia]*lhsc[5,3,b1,b2,b3,b4]
               lhsb[5,3,ic] -= lhsa[5,1,ia]*lhsc[1,3,b1,b2,b3,b4] +
                                 lhsa[5,2,ia]*lhsc[2,3,b1,b2,b3,b4] +
                                 lhsa[5,3,ia]*lhsc[3,3,b1,b2,b3,b4] +
                                 lhsa[5,4,ia]*lhsc[4,3,b1,b2,b3,b4] +
                                 lhsa[5,5,ia]*lhsc[5,3,b1,b2,b3,b4]
               lhsb[1,4,ic] -= lhsa[1,1,ia]*lhsc[1,4,b1,b2,b3,b4] +
                                 lhsa[1,2,ia]*lhsc[2,4,b1,b2,b3,b4] +
                                 lhsa[1,3,ia]*lhsc[3,4,b1,b2,b3,b4] +
                                 lhsa[1,4,ia]*lhsc[4,4,b1,b2,b3,b4] +
                                 lhsa[1,5,ia]*lhsc[5,4,b1,b2,b3,b4]
               lhsb[2,4,ic] -= lhsa[2,1,ia]*lhsc[1,4,b1,b2,b3,b4] +
                                 lhsa[2,2,ia]*lhsc[2,4,b1,b2,b3,b4] +
                                 lhsa[2,3,ia]*lhsc[3,4,b1,b2,b3,b4] +
                                 lhsa[2,4,ia]*lhsc[4,4,b1,b2,b3,b4] +
                                 lhsa[2,5,ia]*lhsc[5,4,b1,b2,b3,b4]
               lhsb[3,4,ic] -= lhsa[3,1,ia]*lhsc[1,4,b1,b2,b3,b4] +
                                 lhsa[3,2,ia]*lhsc[2,4,b1,b2,b3,b4] +
                                 lhsa[3,3,ia]*lhsc[3,4,b1,b2,b3,b4] +
                                 lhsa[3,4,ia]*lhsc[4,4,b1,b2,b3,b4] +
                                 lhsa[3,5,ia]*lhsc[5,4,b1,b2,b3,b4]
               lhsb[4,4,ic] -= lhsa[4,1,ia]*lhsc[1,4,b1,b2,b3,b4] +
                                 lhsa[4,2,ia]*lhsc[2,4,b1,b2,b3,b4] +
                                 lhsa[4,3,ia]*lhsc[3,4,b1,b2,b3,b4] +
                                 lhsa[4,4,ia]*lhsc[4,4,b1,b2,b3,b4] +
                                 lhsa[4,5,ia]*lhsc[5,4,b1,b2,b3,b4]
               lhsb[5,4,ic] -= lhsa[5,1,ia]*lhsc[1,4,b1,b2,b3,b4] +
                                 lhsa[5,2,ia]*lhsc[2,4,b1,b2,b3,b4] +
                                 lhsa[5,3,ia]*lhsc[3,4,b1,b2,b3,b4] +
                                 lhsa[5,4,ia]*lhsc[4,4,b1,b2,b3,b4] +
                                 lhsa[5,5,ia]*lhsc[5,4,b1,b2,b3,b4]
               lhsb[1,5,ic] -= lhsa[1,1,ia]*lhsc[1,5,b1,b2,b3,b4] +
                                 lhsa[1,2,ia]*lhsc[2,5,b1,b2,b3,b4] +
                                 lhsa[1,3,ia]*lhsc[3,5,b1,b2,b3,b4] +
                                 lhsa[1,4,ia]*lhsc[4,5,b1,b2,b3,b4] +
                                 lhsa[1,5,ia]*lhsc[5,5,b1,b2,b3,b4]
               lhsb[2,5,ic] -= lhsa[2,1,ia]*lhsc[1,5,b1,b2,b3,b4] +
                                 lhsa[2,2,ia]*lhsc[2,5,b1,b2,b3,b4] +
                                 lhsa[2,3,ia]*lhsc[3,5,b1,b2,b3,b4] +
                                 lhsa[2,4,ia]*lhsc[4,5,b1,b2,b3,b4] +
                                 lhsa[2,5,ia]*lhsc[5,5,b1,b2,b3,b4]
               lhsb[3,5,ic] -= lhsa[3,1,ia]*lhsc[1,5,b1,b2,b3,b4] +
                                 lhsa[3,2,ia]*lhsc[2,5,b1,b2,b3,b4] +
                                 lhsa[3,3,ia]*lhsc[3,5,b1,b2,b3,b4] +
                                 lhsa[3,4,ia]*lhsc[4,5,b1,b2,b3,b4] +
                                 lhsa[3,5,ia]*lhsc[5,5,b1,b2,b3,b4]
               lhsb[4,5,ic] -= lhsa[4,1,ia]*lhsc[1,5,b1,b2,b3,b4] +
                                 lhsa[4,2,ia]*lhsc[2,5,b1,b2,b3,b4] +
                                 lhsa[4,3,ia]*lhsc[3,5,b1,b2,b3,b4] +
                                 lhsa[4,4,ia]*lhsc[4,5,b1,b2,b3,b4] +
                                 lhsa[4,5,ia]*lhsc[5,5,b1,b2,b3,b4]
               lhsb[5,5,ic] -= lhsa[5,1,ia]*lhsc[1,5,b1,b2,b3,b4] +
                                 lhsa[5,2,ia]*lhsc[2,5,b1,b2,b3,b4] +
                                 lhsa[5,3,ia]*lhsc[3,5,b1,b2,b3,b4] +
                                 lhsa[5,4,ia]*lhsc[4,5,b1,b2,b3,b4] +
                                 lhsa[5,5,ia]*lhsc[5,5,b1,b2,b3,b4] =#


#---------------------------------------------------------------------
#     multiply rhs() by b_inverse() and copy to rhs
#---------------------------------------------------------------------
#               binvrhs(view(lhsb, 1:5, 1:5, isize), view(rhs, 1:5, isize, j, k, c))
               binvrhs(lhsb, isize, rhs, isize, j, k, c)

#=               i = isize
               ir1, ir2, ir3, ir4 = isize, j, k, c

               pivot = 1.00e0/lhsb[1,1,i]
               lhsb[1,2,i] *= pivot
               lhsb[1,3,i] *= pivot
               lhsb[1,4,i] *= pivot
               lhsb[1,5,i] *= pivot
               rhs[1,ir1,ir2,ir3,ir4] *= pivot
         
               coeff = lhsb[2,1,i]
               lhsb[2,2,i] -= coeff*lhsb[1,2,i]
               lhsb[2,3,i] -= coeff*lhsb[1,3,i]
               lhsb[2,4,i] -= coeff*lhsb[1,4,i]
               lhsb[2,5,i] -= coeff*lhsb[1,5,i]
               rhs[2,ir1,ir2,ir3,ir4] -= coeff*rhs[1,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[3,1,i]
               lhsb[3,2,i] -= coeff*lhsb[1,2,i]
               lhsb[3,3,i] -= coeff*lhsb[1,3,i]
               lhsb[3,4,i] -= coeff*lhsb[1,4,i]
               lhsb[3,5,i] -= coeff*lhsb[1,5,i]
               rhs[3,ir1,ir2,ir3,ir4] -= coeff*rhs[1,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[4,1,i]
               lhsb[4,2,i] -= coeff*lhsb[1,2,i]
               lhsb[4,3,i] -= coeff*lhsb[1,3,i]
               lhsb[4,4,i] -= coeff*lhsb[1,4,i]
               lhsb[4,5,i] -= coeff*lhsb[1,5,i]
               rhs[4,ir1,ir2,ir3,ir4] -= coeff*rhs[1,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[5,1,i]
               lhsb[5,2,i] -= coeff*lhsb[1,2,i]
               lhsb[5,3,i] -= coeff*lhsb[1,3,i]
               lhsb[5,4,i] -= coeff*lhsb[1,4,i]
               lhsb[5,5,i] -= coeff*lhsb[1,5,i]
               rhs[5,ir1,ir2,ir3,ir4] -= coeff*rhs[1,ir1,ir2,ir3,ir4]
         
         
               pivot = 1.00e0/lhsb[2,2,i]
               lhsb[2,3,i] *= pivot
               lhsb[2,4,i] *= pivot
               lhsb[2,5,i] *= pivot
               rhs[2,ir1,ir2,ir3,ir4] *= pivot
         
               coeff = lhsb[1,2,i]
               lhsb[1,3,i] -= coeff*lhsb[2,3,i]
               lhsb[1,4,i] -= coeff*lhsb[2,4,i]
               lhsb[1,5,i] -= coeff*lhsb[2,5,i]
               rhs[1,ir1,ir2,ir3,ir4] -= coeff*rhs[2,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[3,2,i]
               lhsb[3,3,i] -= coeff*lhsb[2,3,i]
               lhsb[3,4,i] -= coeff*lhsb[2,4,i]
               lhsb[3,5,i] -= coeff*lhsb[2,5,i]
               rhs[3,ir1,ir2,ir3,ir4] -= coeff*rhs[2,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[4,2,i]
               lhsb[4,3,i] -= coeff*lhsb[2,3,i]
               lhsb[4,4,i] -= coeff*lhsb[2,4,i]
               lhsb[4,5,i] -= coeff*lhsb[2,5,i]
               rhs[4,ir1,ir2,ir3,ir4] -= coeff*rhs[2,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[5,2,i]
               lhsb[5,3,i] -= coeff*lhsb[2,3,i]
               lhsb[5,4,i] -= coeff*lhsb[2,4,i]
               lhsb[5,5,i] -= coeff*lhsb[2,5,i]
               rhs[5,ir1,ir2,ir3,ir4] -= coeff*rhs[2,ir1,ir2,ir3,ir4]
         
         
               pivot = 1.00e0/lhsb[3,3,i]
               lhsb[3,4,i] *= pivot
               lhsb[3,5,i] *= pivot
               rhs[3,ir1,ir2,ir3,ir4] *= pivot
         
               coeff = lhsb[1,3,i]
               lhsb[1,4,i] -= coeff*lhsb[3,4,i]
               lhsb[1,5,i] -= coeff*lhsb[3,5,i]
               rhs[1,ir1,ir2,ir3,ir4] -= coeff*rhs[3,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[2,3,i]
               lhsb[2,4,i] -= coeff*lhsb[3,4,i]
               lhsb[2,5,i] -= coeff*lhsb[3,5,i]
               rhs[2,ir1,ir2,ir3,ir4]  -= coeff*rhs[3,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[4,3,i]
               lhsb[4,4,i] -= coeff*lhsb[3,4,i]
               lhsb[4,5,i] -= coeff*lhsb[3,5,i]
               rhs[4,ir1,ir2,ir3,ir4] -= coeff*rhs[3,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[5,3,i]
               lhsb[5,4,i] -= coeff*lhsb[3,4,i]
               lhsb[5,5,i] -= coeff*lhsb[3,5,i]
               rhs[5,ir1,ir2,ir3,ir4] -= coeff*rhs[3,ir1,ir2,ir3,ir4]
         
         
               pivot = 1.00e0/lhsb[4,4,i]
               lhsb[4,5,i] *= pivot
               rhs[4,ir1,ir2,ir3,ir4] *= pivot
         
               coeff = lhsb[1,4,i]
               lhsb[1,5,i] -= coeff*lhsb[4,5,i]
               rhs[1,ir1,ir2,ir3,ir4] -= coeff*rhs[4,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[2,4,i]
               lhsb[2,5,i] -= coeff*lhsb[4,5,i]
               rhs[2,ir1,ir2,ir3,ir4] -= coeff*rhs[4,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[3,4,i]
               lhsb[3,5,i] -= coeff*lhsb[4,5,i]
               rhs[3,ir1,ir2,ir3,ir4] -= coeff*rhs[4,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[5,4,i]
               lhsb[5,5,i] -= coeff*lhsb[4,5,i]
               rhs[5,ir1,ir2,ir3,ir4] -= coeff*rhs[4,ir1,ir2,ir3,ir4]
               
         
               pivot = 1.00e0/lhsb[5,5,i]
               rhs[5,ir1,ir2,ir3,ir4] *= pivot
         
               coeff = lhsb[1,5,i]
               rhs[1,ir1,ir2,ir3,ir4] -= coeff*rhs[5,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[2,5,i]
               rhs[2,ir1,ir2,ir3,ir4] -= coeff*rhs[5,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[3,5,i]
               rhs[3,ir1,ir2,ir3,ir4] -= coeff*rhs[5,ir1,ir2,ir3,ir4]
         
               coeff = lhsb[4,5,i]
               rhs[4,ir1,ir2,ir3,ir4] -= coeff*rhs[5,ir1,ir2,ir3,ir4] =#
         
            end
         end
      end


      return nothing
   end
end

