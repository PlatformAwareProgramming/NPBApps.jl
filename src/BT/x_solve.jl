using FortranFiles
using OffsetArrays
using Parameters
using Printf


#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function x_solve()

#---------------------------------------------------------------------
#---------------------------------------------------------------------

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

#      use bt_data
#      use mpinpb

#      implicit none

#      integer  c, istart, stage,  
#           FIRST, LAST, recv_id, ERROR, r_status[MPI_STATUS_SIZE],  
#           isize,jsize,ksize,send_id

      istart = 0

      if (timeron) timer_start(t_xsolve) end
#---------------------------------------------------------------------
#     in our terminology stage is the number of the cell in the x-direction
#     i.e. stage = 1 means the start of the line stage=ncells means end
#---------------------------------------------------------------------
      for stage = 1:ncells
         c = slice[1, stage]
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
#            call lhsx(c)
            x_solve_cell(FIRST, LAST, c)
         else
#---------------------------------------------------------------------
#     Not the first cell of this line, so receive info from
#     processor working on preceeding cell
#---------------------------------------------------------------------
            FIRST = 0
            if (timeron) timer_start(t_xcomm) end
            x_receive_solve_info(recv_id, c)
#---------------------------------------------------------------------
#     overlap computations and communications
#---------------------------------------------------------------------
#            call lhsx(c)
#---------------------------------------------------------------------
#     wait for completion
#---------------------------------------------------------------------
            mpi_wait(send_id, r_status, ERROR)
            mpi_wait(recv_id, r_status, ERROR)
            if (timeron) timer_stop(t_xcomm) end
#---------------------------------------------------------------------
#     install C'(istart) and rhs'(istart) to be used in this cell
#---------------------------------------------------------------------
            x_unpack_solve_info(c)
            x_solve_cell(FIRST, LAST, c)
         end

         if (LAST == 0) x_send_solve_info(send_id, c) end
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
            x_backsubstitute(FIRST, LAST, c)
         else
            if (timeron) timer_start(t_xcomm) end
            x_receive_backsub_info(recv_id, c)
            mpi_wait(send_id, r_status, ERROR)
            mpi_wait(recv_id, r_status, ERROR)
            if (timeron) timer_stop(t_xcomm) end
            x_unpack_backsub_info(c)
            x_backsubstitute(FIRST, LAST, c)
         end
         if (FIRST == 0) x_send_backsub_info(send_id, c) end
      end

      if (timeron) timer_stop(t_xsolve) end

      return nothing
      end


#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function x_unpack_solve_info(c)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     unpack C'(-1) and rhs'(-1) for
#     all j and k
#---------------------------------------------------------------------

#      use bt_data
#      implicit none

#      integer j,k,m,n,ptr,c,istart

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

#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function x_send_solve_info(send_id, c)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     pack up and send C'(iend) and rhs'(iend) for
#     all j and k
#---------------------------------------------------------------------

#      use bt_data
#      use mpinpb

#      implicit none

#      integer j,k,m,n,isize,ptr,c,jp,kp
#      integer ERROR,send_id,buffer_size

      isize = cell_size[1, c]-1
      jp = cell_coord[2, c] - 1
      kp = cell_coord[3, c] - 1
      buffer_size = MAX_CELL_DIM*MAX_CELL_DIM*(
           BLOCK_SIZE*BLOCK_SIZE + BLOCK_SIZE)

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
      if (timeron) timer_start(t_xcomm) end
      MPI.Isend(in_buffer, buffer_size,
           dp_type, successor[1],
           WEST+jp+kp*NCELLS, comm_solve,
           send_id, ERROR)
      if (timeron) timer_stop(t_xcomm) end

      return nothing
      end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function x_send_backsub_info(send_id, c)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     pack up and send u[istart] for all j and k
#---------------------------------------------------------------------

#      use bt_data
#      use mpinpb

#      implicit none

#      integer j,k,n,ptr,c,istart,jp,kp
#      integer ERROR,send_id,buffer_size

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
      if (timeron) timer_start(t_xcomm) end
      MPI.Isend(in_buffer, buffer_size,
           dp_type, predecessor[1],
           EAST+jp+kp*NCELLS, comm_solve,
           send_id, ERROR)
      if (timeron) timer_stop(t_xcomm) end

      return nothing
      end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function x_unpack_backsub_info(c)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     unpack u[isize] for all j and k
#---------------------------------------------------------------------

#      use bt_data
#      implicit none

#      integer j,k,n,ptr,c

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

#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function x_receive_backsub_info(recv_id, c)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     post mpi receives
#---------------------------------------------------------------------

#      use bt_data
#      use mpinpb

#      implicit none

#      integer ERROR,recv_id,jp,kp,c,buffer_size

      jp = cell_coord[2, c] - 1
      kp = cell_coord[3, c] - 1
      buffer_size = MAX_CELL_DIM*MAX_CELL_DIM*BLOCK_SIZE
      MPI.Irecv!(out_buffer, buffer_size,
           dp_type, successor[1],
           EAST+jp+kp*NCELLS, comm_solve,
           recv_id, ERROR)

      return nothing
      end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function x_receive_solve_info(recv_id, c)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     post mpi receives 
#---------------------------------------------------------------------

#      use bt_data
#      use mpinpb

#      implicit none

#      integer jp,kp,recv_id,ERROR,c,buffer_size

      jp = cell_coord[2, c] - 1
      kp = cell_coord[3, c] - 1
      buffer_size = MAX_CELL_DIM*MAX_CELL_DIM*(
           BLOCK_SIZE*BLOCK_SIZE + BLOCK_SIZE)
      MPI.Irecv!(out_buffer, buffer_size,
           dp_type, predecessor[1],
           WEST+jp+kp*NCELLS,  comm_solve,
           recv_id, ERROR)

      return nothing
      end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function x_backsubstitute(FIRST, LAST, c)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     back solve: if last cell, then generate u[isize]=rhs[isize]
#     else assume u[isize] is loaded in un pack backsub_info
#     so just use it
#     after call u[istart] will be sent to next cell
#---------------------------------------------------------------------

#      use bt_data
#      implicit none

#      integer FIRST, LAST, c, i, j, k
#      integer m,n,isize,jsize,ksize,istart

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
                     rhs[m, isize, j, k, c] = rhs[m, isize, j, k, c]-
                            lhsc[m, n, isize, j, k, c]*
                          backsub_info[n, j, k, c]
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
                     rhs[m, i, j, k, c] = rhs[m, i, j, k, c]-
                            lhsc[m, n, i, j, k, c]*rhs[n, i+1, j, k, c]
                  end
               end
            end
         end
      end

      return nothing
      end


#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function x_solve_cell(FIRST, LAST, c)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     performs guaussian elimination on this cell.
#     
#     assumes that unpacking routines for non-first cells 
#     preload C' and rhs' from previous cell.
#     
#     assumed send happens outside this routine, but that
#     c'(IMAX) and rhs'(IMAX) will be sent to next cell
#---------------------------------------------------------------------

#      use bt_data
#      implicit none

#      DOUBLEPRECISION tmp1, tmp2, tmp3
#      integer FIRST,LAST,c
#      integer i,j,k,isize,ksize,jsize,istart

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

               tmp1 = rho_i(i, j, k, c)
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

               fjac[2, 1, i] = -(u[2, i, j, k, c] * tmp2 *
                    u[2, i, j, k, c])+
                     c2 * qs[i, j, k, c]
               fjac[2, 2, i] = ( 2.0e+00 - c2 )*
                     ( u[2, i, j, k, c] * tmp1 )
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

               fjac[5, 1, i] = ( c2 * 2.0e0 * qs[i, j, k, c]-
                     c1 * ( u[5, i, j, k, c] * tmp1 ) )*
                     ( u[2, i, j, k, c] * tmp1 )
               fjac[5, 2, i] = c1 *  u[5, i, j, k, c] * tmp1-
                     c2*
                     ( u[2, i, j, k, c]*u[2, i, j, k, c] * tmp2+
                     qs[i, j, k, c] )
               fjac[5, 3, i] = - c2 * ( u[3, i, j, k, c]*u[2, i, j, k, c] )*
                     tmp2
               fjac[5, 4, i] = - c2 * ( u[4, i, j, k, c]*u[2, i, j, k, c] )*
                     tmp2
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

               njac[5, 2, i] = ( con43 * c3c4-
                     c1345 ) * tmp2 * u[2, i, j, k, c]
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

                lhsa[1, 1, i] = - tmp2 * fjac[1, 1, i-1]-
                     tmp1 * njac[1, 1, i-1]-
                     tmp1 * dx1
                lhsa[1, 2, i] = - tmp2 * fjac[1, 2, i-1]-
                     tmp1 * njac[1, 2, i-1]
                lhsa[1, 3, i] = - tmp2 * fjac[1, 3, i-1]-
                     tmp1 * njac[1, 3, i-1]
                lhsa[1, 4, i] = - tmp2 * fjac[1, 4, i-1]-
                     tmp1 * njac[1, 4, i-1]
                lhsa[1, 5, i] = - tmp2 * fjac[1, 5, i-1]-
                     tmp1 * njac[1, 5, i-1]

                lhsa[2, 1, i] = - tmp2 * fjac[2, 1, i-1]-
                     tmp1 * njac[2, 1, i-1]
                lhsa[2, 2, i] = - tmp2 * fjac[2, 2, i-1]-
                     tmp1 * njac[2, 2, i-1]-
                     tmp1 * dx2
                lhsa[2, 3, i] = - tmp2 * fjac[2, 3, i-1]-
                     tmp1 * njac[2, 3, i-1]
                lhsa[2, 4, i] = - tmp2 * fjac[2, 4, i-1]-
                     tmp1 * njac[2, 4, i-1]
                lhsa[2, 5, i] = - tmp2 * fjac[2, 5, i-1]-
                     tmp1 * njac[2, 5, i-1]

                lhsa[3, 1, i] = - tmp2 * fjac[3, 1, i-1]-
                     tmp1 * njac[3, 1, i-1]
                lhsa[3, 2, i] = - tmp2 * fjac[3, 2, i-1]-
                     tmp1 * njac[3, 2, i-1]
                lhsa[3, 3, i] = - tmp2 * fjac[3, 3, i-1]-
                     tmp1 * njac[3, 3, i-1]-
                     tmp1 * dx3
                lhsa[3, 4, i] = - tmp2 * fjac[3, 4, i-1]-
                     tmp1 * njac[3, 4, i-1]
                lhsa[3, 5, i] = - tmp2 * fjac[3, 5, i-1]-
                     tmp1 * njac[3, 5, i-1]

                lhsa[4, 1, i] = - tmp2 * fjac[4, 1, i-1]-
                     tmp1 * njac[4, 1, i-1]
                lhsa[4, 2, i] = - tmp2 * fjac[4, 2, i-1]-
                     tmp1 * njac[4, 2, i-1]
                lhsa[4, 3, i] = - tmp2 * fjac[4, 3, i-1]-
                     tmp1 * njac[4, 3, i-1]
                lhsa[4, 4, i] = - tmp2 * fjac[4, 4, i-1]-
                     tmp1 * njac[4, 4, i-1]-
                     tmp1 * dx4
                lhsa[4, 5, i] = - tmp2 * fjac[4, 5, i-1]-
                     tmp1 * njac[4, 5, i-1]

                lhsa[5, 1, i] = - tmp2 * fjac[5, 1, i-1]-
                     tmp1 * njac[5, 1, i-1]
                lhsa[5, 2, i] = - tmp2 * fjac[5, 2, i-1]-
                     tmp1 * njac[5, 2, i-1]
                lhsa[5, 3, i] = - tmp2 * fjac[5, 3, i-1]-
                     tmp1 * njac[5, 3, i-1]
                lhsa[5, 4, i] = - tmp2 * fjac[5, 4, i-1]-
                     tmp1 * njac[5, 4, i-1]
                lhsa[5, 5, i] = - tmp2 * fjac[5, 5, i-1]-
                     tmp1 * njac[5, 5, i-1]-
                     tmp1 * dx5

                lhsb[1, 1, i] = 1.0e+00+
                     tmp1 * 2.0e+00 * njac[1, 1, i]+
                     tmp1 * 2.0e+00 * dx1
                lhsb[1, 2, i] = tmp1 * 2.0e+00 * njac[1, 2, i]
                lhsb[1, 3, i] = tmp1 * 2.0e+00 * njac[1, 3, i]
                lhsb[1, 4, i] = tmp1 * 2.0e+00 * njac[1, 4, i]
                lhsb[1, 5, i] = tmp1 * 2.0e+00 * njac[1, 5, i]

                lhsb[2, 1, i] = tmp1 * 2.0e+00 * njac[2, 1, i]
                lhsb[2, 2, i] = 1.0e+00+
                     tmp1 * 2.0e+00 * njac[2, 2, i]+
                     tmp1 * 2.0e+00 * dx2
                lhsb[2, 3, i] = tmp1 * 2.0e+00 * njac[2, 3, i]
                lhsb[2, 4, i] = tmp1 * 2.0e+00 * njac[2, 4, i]
                lhsb[2, 5, i] = tmp1 * 2.0e+00 * njac[2, 5, i]

                lhsb[3, 1, i] = tmp1 * 2.0e+00 * njac[3, 1, i]
                lhsb[3, 2, i] = tmp1 * 2.0e+00 * njac[3, 2, i]
                lhsb[3, 3, i] = 1.0e+00+
                     tmp1 * 2.0e+00 * njac[3, 3, i]+
                     tmp1 * 2.0e+00 * dx3
                lhsb[3, 4, i] = tmp1 * 2.0e+00 * njac[3, 4, i]
                lhsb[3, 5, i] = tmp1 * 2.0e+00 * njac[3, 5, i]

                lhsb[4, 1, i] = tmp1 * 2.0e+00 * njac[4, 1, i]
                lhsb[4, 2, i] = tmp1 * 2.0e+00 * njac[4, 2, i]
                lhsb[4, 3, i] = tmp1 * 2.0e+00 * njac[4, 3, i]
                lhsb[4, 4, i] = 1.0e+00+
                     tmp1 * 2.0e+00 * njac[4, 4, i]+
                     tmp1 * 2.0e+00 * dx4
                lhsb[4, 5, i] = tmp1 * 2.0e+00 * njac[4, 5, i]

                lhsb[5, 1, i] = tmp1 * 2.0e+00 * njac[5, 1, i]
                lhsb[5, 2, i] = tmp1 * 2.0e+00 * njac[5, 2, i]
                lhsb[5, 3, i] = tmp1 * 2.0e+00 * njac[5, 3, i]
                lhsb[5, 4, i] = tmp1 * 2.0e+00 * njac[5, 4, i]
                lhsb[5, 5, i] = 1.0e+00+
                     tmp1 * 2.0e+00 * njac[5, 5, i]+
                     tmp1 * 2.0e+00 * dx5

                lhsc[1, 1, i, j, k, c] =  tmp2 * fjac[1, 1, i+1]-
                     tmp1 * njac[1, 1, i+1]-
                     tmp1 * dx1
                lhsc[1, 2, i, j, k, c] =  tmp2 * fjac[1, 2, i+1]-
                     tmp1 * njac[1, 2, i+1]
                lhsc[1, 3, i, j, k, c] =  tmp2 * fjac[1, 3, i+1]-
                     tmp1 * njac[1, 3, i+1]
                lhsc[1, 4, i, j, k, c] =  tmp2 * fjac[1, 4, i+1]-
                     tmp1 * njac[1, 4, i+1]
                lhsc[1, 5, i, j, k, c] =  tmp2 * fjac[1, 5, i+1]-
                     tmp1 * njac[1, 5, i+1]

                lhsc[2, 1, i, j, k, c] =  tmp2 * fjac[2, 1, i+1]-
                     tmp1 * njac[2, 1, i+1]
                lhsc[2, 2, i, j, k, c] =  tmp2 * fjac[2, 2, i+1]-
                     tmp1 * njac[2, 2, i+1]-
                     tmp1 * dx2
                lhsc[2, 3, i, j, k, c] =  tmp2 * fjac[2, 3, i+1]-
                     tmp1 * njac[2, 3, i+1]
                lhsc[2, 4, i, j, k, c] =  tmp2 * fjac[2, 4, i+1]-
                     tmp1 * njac[2, 4, i+1]
                lhsc[2, 5, i, j, k, c] =  tmp2 * fjac[2, 5, i+1]-
                     tmp1 * njac[2, 5, i+1]

                lhsc[3, 1, i, j, k, c] =  tmp2 * fjac[3, 1, i+1]-
                     tmp1 * njac[3, 1, i+1]
                lhsc[3, 2, i, j, k, c] =  tmp2 * fjac[3, 2, i+1]-
                     tmp1 * njac[3, 2, i+1]
                lhsc[3, 3, i, j, k, c] =  tmp2 * fjac[3, 3, i+1]-
                     tmp1 * njac[3, 3, i+1]-
                     tmp1 * dx3
                lhsc[3, 4, i, j, k, c] =  tmp2 * fjac[3, 4, i+1]-
                     tmp1 * njac[3, 4, i+1]
                lhsc[3, 5, i, j, k, c] =  tmp2 * fjac[3, 5, i+1]-
                     tmp1 * njac[3, 5, i+1]

                lhsc[4, 1, i, j, k, c] =  tmp2 * fjac[4, 1, i+1]-
                     tmp1 * njac[4, 1, i+1]
                lhsc[4, 2, i, j, k, c] =  tmp2 * fjac[4, 2, i+1]-
                     tmp1 * njac[4, 2, i+1]
                lhsc[4, 3, i, j, k, c] =  tmp2 * fjac[4, 3, i+1]-
                     tmp1 * njac[4, 3, i+1]
                lhsc[4, 4, i, j, k, c] =  tmp2 * fjac[4, 4, i+1]-
                     tmp1 * njac[4, 4, i+1]-
                     tmp1 * dx4
                lhsc[4, 5, i, j, k, c] =  tmp2 * fjac[4, 5, i+1]-
                     tmp1 * njac[4, 5, i+1]

                lhsc[5, 1, i, j, k, c] =  tmp2 * fjac[5, 1, i+1]-
                     tmp1 * njac[5, 1, i+1]
                lhsc[5, 2, i, j, k, c] =  tmp2 * fjac[5, 2, i+1]-
                     tmp1 * njac[5, 2, i+1]
                lhsc[5, 3, i, j, k, c] =  tmp2 * fjac[5, 3, i+1]-
                     tmp1 * njac[5, 3, i+1]
                lhsc[5, 4, i, j, k, c] =  tmp2 * fjac[5, 4, i+1]-
                     tmp1 * njac[5, 4, i+1]
                lhsc[5, 5, i, j, k, c] =  tmp2 * fjac[5, 5, i+1]-
                     tmp1 * njac[5, 5, i+1]-
                     tmp1 * dx5

            end


#---------------------------------------------------------------------
#     outer most do loops - sweeping in i direction
#---------------------------------------------------------------------
            if FIRST == 1

#---------------------------------------------------------------------
#     multiply c[istart,j,k] by b_inverse and copy back to c
#     multiply rhs[istart] by b_inverse(istart) and copy to rhs
#---------------------------------------------------------------------
               binvcrhs(  lhsb[1, 1, istart],
                               lhsc[1, 1, istart, j, k, c],
                              rhs[1, istart, j, k, c] )

            end

#---------------------------------------------------------------------
#     begin inner most do loop
#     do all the elements of the cell unless last 
#---------------------------------------------------------------------
            for i = istart+FIRST:isize-LAST

#---------------------------------------------------------------------
#     rhs[i] = rhs[i] - A*rhs[i-1]
#---------------------------------------------------------------------
               matvec_sub( lhsa[1, 1, i],
                               rhs[1, i-1, j, k, c], rhs[1, i, j, k, c])

#---------------------------------------------------------------------
#      b[i] =  b[i] - C(i-1)*A(i)
#---------------------------------------------------------------------
               matmul_sub( lhsa[1, 1, i],
                                lhsc[1, 1, i-1, j, k, c],
                                lhsb[1, 1, i])


#---------------------------------------------------------------------
#     multiply c[i,j,k] by b_inverse and copy back to c
#     multiply rhs[1,j,k] by b_inverse(1,j,k) and copy to rhs
#---------------------------------------------------------------------
               binvcrhs(  lhsb[1, 1, i],
                               lhsc[1, 1, i, j, k, c],
                              rhs[1, i, j, k, c] )

            end

#---------------------------------------------------------------------
#     Now finish up special cases for last cell
#---------------------------------------------------------------------
            if LAST == 1

#---------------------------------------------------------------------
#     rhs[isize] = rhs[isize] - A*rhs[isize-1]
#---------------------------------------------------------------------
               matvec_sub( lhsa[1, 1, isize],
                               rhs[1, isize-1, j, k, c], rhs[1, isize, j, k, c])

#---------------------------------------------------------------------
#      b[isize] =  b[isize] - C(isize-1)*A(isize)
#---------------------------------------------------------------------
               matmul_sub( lhsa[1, 1, isize],
                                lhsc[1, 1, isize-1, j, k, c],
                                lhsb[1, 1, isize])

#---------------------------------------------------------------------
#     multiply rhs() by b_inverse() and copy to rhs
#---------------------------------------------------------------------
               binvrhs(  lhsb[1, 1, isize],
                             rhs[1, isize, j, k, c] )

            end
         end
      end


      return nothing
      end

