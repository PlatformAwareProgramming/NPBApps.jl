using FortranFiles
using OffsetArrays
using Parameters
using Printf

#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function z_solve()

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     Performs line solves in Z direction by first factoring
#     the block-tridiagonal matrix into an upper triangular matrix, 
#     and then performing back substitution to solve for the unknow
#     vectors of each line.  
#     
#     Make sure we treat elements zero to cell_size in the direction
#     of the sweep.
#---------------------------------------------------------------------

#      use bt_data
#      use mpinpb

#      implicit none

#      integer c, kstart, stage,  
#           FIRST, LAST, recv_id, ERROR, r_status[MPI_STATUS_SIZE],  
#           isize,jsize,ksize,send_id

      kstart = 0

      if (timeron) timer_start(t_zsolve) end
#---------------------------------------------------------------------
#     in our terminology stage is the number of the cell in the y-direct
#     i.e. stage = 1 means the start of the line stage=ncells means end
#---------------------------------------------------------------------
      for stage = 1:ncells
         c = slice[3, stage]
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
#            call lhsz(c)
            z_solve_cell(FIRST, LAST, c)
         else
#---------------------------------------------------------------------
#     Not the first cell of this line, so receive info from
#     processor working on preceeding cell
#---------------------------------------------------------------------
            FIRST = 0
            if (timeron) timer_start(t_zcomm) end
            z_receive_solve_info(recv_id, c)
#---------------------------------------------------------------------
#     overlap computations and communications
#---------------------------------------------------------------------
#            call lhsz(c)
#---------------------------------------------------------------------
#     wait for completion
#---------------------------------------------------------------------
            mpi_wait(send_id, r_status, ERROR)
            mpi_wait(recv_id, r_status, ERROR)
            if (timeron) timer_stop(t_zcomm) end
#---------------------------------------------------------------------
#     install C'(kstart+1) and rhs'(kstart+1) to be used in this cell
#---------------------------------------------------------------------
            z_unpack_solve_info(c)
            z_solve_cell(FIRST, LAST, c)
         end

         if (LAST == 0) z_send_solve_info(send_id, c) end
      end

#---------------------------------------------------------------------
#     now perform backsubstitution in reverse direction
#---------------------------------------------------------------------
      for stage = ncells:-1:1
         c = slice[3, stage]
         FIRST = 0
         LAST = 0
         if (stage == 1) FIRST = 1 end
         if stage == ncells
            LAST = 1
#---------------------------------------------------------------------
#     last cell, so perform back substitute without waiting
#---------------------------------------------------------------------
            z_backsubstitute(FIRST, LAST, c)
         else
            if (timeron) timer_start(t_zcomm) end
            z_receive_backsub_info(recv_id, c)
            mpi_wait(send_id, r_status, ERROR)
            mpi_wait(recv_id, r_status, ERROR)
            if (timeron) timer_stop(t_zcomm) end
            z_unpack_backsub_info(c)
            z_backsubstitute(FIRST, LAST, c)
         end
         if (FIRST == 0) z_send_backsub_info(send_id, c) end
      end

      if (timeron) timer_stop(t_zsolve) end

      return nothing
      end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function z_unpack_solve_info(c)
#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     unpack C'(-1) and rhs'(-1) for
#     all i and j
#---------------------------------------------------------------------

#      use bt_data
#      implicit none

#      integer i,j,m,n,ptr,c,kstart

      kstart = 0
      ptr = 0
      for j = 0:JMAX-1
         for i = 0:IMAX-1
            for m = 1:BLOCK_SIZE
               for n = 1:BLOCK_SIZE
                   lhsc[m, n, i, j, kstart-1, c] = out_buffer[ptr+n]
               end
               ptr = ptr+BLOCK_SIZE
            end
            for n = 1:BLOCK_SIZE
               rhs[n, i, j, kstart-1, c] = out_buffer[ptr+n]
            end
            ptr = ptr+BLOCK_SIZE
         end
      end

      return nothing
      end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function z_send_solve_info(send_id, c)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     pack up and send C'(kend) and rhs'(kend) for
#     all i and j
#---------------------------------------------------------------------

#      use bt_data
#      use mpinpb

#      implicit none

#      integer i,j,m,n,ksize,ptr,c,ip,jp
#      integer ERROR,send_id,buffer_size

      ksize = cell_size[3, c]-1
      ip = cell_coord[1, c] - 1
      jp = cell_coord[2, c] - 1
      buffer_size = MAX_CELL_DIM*MAX_CELL_DIM*(
           BLOCK_SIZE*BLOCK_SIZE + BLOCK_SIZE)

#---------------------------------------------------------------------
#     pack up buffer
#---------------------------------------------------------------------
      ptr = 0
      for j = 0:JMAX-1
         for i = 0:IMAX-1
            for m = 1:BLOCK_SIZE
               for n = 1:BLOCK_SIZE
                  in_buffer[ptr+n] =  lhsc[m, n, i, j, ksize, c]
               end
               ptr = ptr+BLOCK_SIZE
            end
            for n = 1:BLOCK_SIZE
               in_buffer[ptr+n] = rhs[n, i, j, ksize, c]
            end
            ptr = ptr+BLOCK_SIZE
         end
      end

#---------------------------------------------------------------------
#     send buffer 
#---------------------------------------------------------------------
      if (timeron) timer_start(t_zcomm) end
      MPI.Isend(in_buffer, buffer_size,
           dp_type, successor[3],
           BOTTOM+ip+jp*NCELLS, comm_solve,
           send_id, ERROR)
      if (timeron) timer_stop(t_zcomm) end

      return nothing
      end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function z_send_backsub_info(send_id, c)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     pack up and send u[jstart] for all i and j
#---------------------------------------------------------------------

#      use bt_data
#      use mpinpb

#      implicit none

#      integer i,j,n,ptr,c,kstart,ip,jp
#      integer ERROR,send_id,buffer_size

#---------------------------------------------------------------------
#     Send element 0 to previous processor
#---------------------------------------------------------------------
      kstart = 0
      ip = cell_coord[1, c]-1
      jp = cell_coord[2, c]-1
      buffer_size = MAX_CELL_DIM*MAX_CELL_DIM*BLOCK_SIZE
      ptr = 0
      for j = 0:JMAX-1
         for i = 0:IMAX-1
            for n = 1:BLOCK_SIZE
               in_buffer[ptr+n] = rhs[n, i, j, kstart, c]
            end
            ptr = ptr+BLOCK_SIZE
         end
      end

      if (timeron) timer_start(t_zcomm) end
      MPI.Isend(in_buffer, buffer_size,
           dp_type, predecessor[3],
           TOP+ip+jp*NCELLS, comm_solve,
           send_id, ERROR)
      if (timeron) timer_stop(t_zcomm) end

      return nothing
      end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function z_unpack_backsub_info(c)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     unpack u[ksize] for all i and j
#---------------------------------------------------------------------

#      use bt_data
#      implicit none

#      integer i,j,n,ptr,c

      ptr = 0
      for j = 0:JMAX-1
         for i = 0:IMAX-1
            for n = 1:BLOCK_SIZE
               backsub_info[n, i, j, c] = out_buffer[ptr+n]
            end
            ptr = ptr+BLOCK_SIZE
         end
      end

      return nothing
      end


#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function z_receive_backsub_info(recv_id, c)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     post mpi receives
#---------------------------------------------------------------------

#      use bt_data
#      use mpinpb

#      implicit none

#      integer ERROR,recv_id,ip,jp,c,buffer_size
      ip = cell_coord[1, c] - 1
      jp = cell_coord[2, c] - 1
      buffer_size = MAX_CELL_DIM*MAX_CELL_DIM*BLOCK_SIZE
      MPI.Irecv!(out_buffer, buffer_size,
           dp_type, successor[3],
           TOP+ip+jp*NCELLS, comm_solve,
           recv_id, ERROR)

      return nothing
      end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function z_receive_solve_info(recv_id, c)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     post mpi receives 
#---------------------------------------------------------------------

#      use bt_data
#      use mpinpb

#      implicit none

#      integer ip,jp,recv_id,ERROR,c,buffer_size
      ip = cell_coord[1, c] - 1
      jp = cell_coord[2, c] - 1
      buffer_size = MAX_CELL_DIM*MAX_CELL_DIM*(
           BLOCK_SIZE*BLOCK_SIZE + BLOCK_SIZE)
      MPI.Irecv!(out_buffer, buffer_size,
           dp_type, predecessor[3],
           BOTTOM+ip+jp*NCELLS, comm_solve,
           recv_id, ERROR)

      return nothing
      end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function z_backsubstitute(FIRST, LAST, c)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     back solve: if last cell, then generate u[ksize]=rhs[ksize]
#     else assume u[ksize] is loaded in un pack backsub_info
#     so just use it
#     after call u[kstart] will be sent to next cell
#---------------------------------------------------------------------

#      use bt_data
#      implicit none

#      integer FIRST, LAST, c, i, k
#      integer m,n,j,jsize,isize,ksize,kstart

      kstart = 0
      isize = cell_size[1, c]-cell_end[1, c]-1
      jsize = cell_size[2, c]-cell_end[2, c]-1
      ksize = cell_size[3, c]-1
      if LAST == 0
         for j = cell_start[2, c]:jsize
            for i = cell_start[1, c]:isize
#---------------------------------------------------------------------
#     u[jsize] uses info from previous cell if not last cell
#---------------------------------------------------------------------
               for m = 1:BLOCK_SIZE
                  for n = 1:BLOCK_SIZE
                     rhs[m, i, j, ksize, c] = rhs[m, i, j, ksize, c]-
                            lhsc[m, n, i, j, ksize, c]*
                          backsub_info[n, i, j, c]
                  end
               end
            end
         end
      end
      for k = ksize-1:-1:kstart
         for j = cell_start[2, c]:jsize
            for i = cell_start[1, c]:isize
               for m = 1:BLOCK_SIZE
                  for n = 1:BLOCK_SIZE
                     rhs[m, i, j, k, c] = rhs[m, i, j, k, c]-
                            lhsc[m, n, i, j, k, c]*rhs[n, i, j, k+1, c]
                  end
               end
            end
         end
      end

      return nothing
      end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function z_solve_cell(FIRST, LAST, c)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     performs guaussian elimination on this cell.
#     
#     assumes that unpacking routines for non-first cells 
#     preload C' and rhs' from previous cell.
#     
#     assumed send happens outside this routine, but that
#     c'(KMAX) and rhs'(KMAX) will be sent to next cell.
#---------------------------------------------------------------------

#      use bt_data
#      implicit none

#      DOUBLEPRECISION tmp1, tmp2, tmp3
#      integer FIRST,LAST,c
#      integer i,j,k,m,n,isize,ksize,jsize,kstart

      kstart = 0
      isize = cell_size[1, c]-cell_end[1, c]-1
      jsize = cell_size[2, c]-cell_end[2, c]-1
      ksize = cell_size[3, c]-1

#---------------------------------------------------------------------
#     zero the left hand side for starters
#     set diagonal values to 1. This is overkill, but convenient
#---------------------------------------------------------------------
      for i = 0:isize
         for m = 1:5
            for n = 1:5
                lhsa[m, n, i, 0] = 0.0e0
                lhsb[m, n, i, 0] = 0.0e0
                lhsa[m, n, i, ksize] = 0.0e0
                lhsb[m, n, i, ksize] = 0.0e0
            end
             lhsb[m, m, i, 0] = 1.0e0
             lhsb[m, m, i, ksize] = 1.0e0
         end
      end

      for j = cell_start[2, c]:jsize

#---------------------------------------------------------------------
#     This function computes the left hand side for the three z-factors 
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     Compute the indices for storing the block-diagonal matrix;
#     determine c (labeled f) and s jacobians for cell !
#---------------------------------------------------------------------

         for k = cell_start[3, c]-1:cell_size[3, c]-cell_end[3, c]
            for i = cell_start[1, c]:isize

               tmp1 = 1.0e0 / u[1, i, j, k, c]
               tmp2 = tmp1 * tmp1
               tmp3 = tmp1 * tmp2

               fjac[1, 1, i, k] = 0.0e+00
               fjac[1, 2, i, k] = 0.0e+00
               fjac[1, 3, i, k] = 0.0e+00
               fjac[1, 4, i, k] = 1.0e+00
               fjac[1, 5, i, k] = 0.0e+00

               fjac[2, 1, i, k] = - ( u[2, i, j, k, c]*u[4, i, j, k, c] )*
                     tmp2
               fjac[2, 2, i, k] = u[4, i, j, k, c] * tmp1
               fjac[2, 3, i, k] = 0.0e+00
               fjac[2, 4, i, k] = u[2, i, j, k, c] * tmp1
               fjac[2, 5, i, k] = 0.0e+00

               fjac[3, 1, i, k] = - ( u[3, i, j, k, c]*u[4, i, j, k, c] )*
                     tmp2
               fjac[3, 2, i, k] = 0.0e+00
               fjac[3, 3, i, k] = u[4, i, j, k, c] * tmp1
               fjac[3, 4, i, k] = u[3, i, j, k, c] * tmp1
               fjac[3, 5, i, k] = 0.0e+00

               fjac[4, 1, i, k] = - (u[4, i, j, k, c]*u[4, i, j, k, c] * tmp2 )+
                     c2 * qs[i, j, k, c]
               fjac[4, 2, i, k] = - c2 *  u[2, i, j, k, c] * tmp1
               fjac[4, 3, i, k] = - c2 *  u[3, i, j, k, c] * tmp1
               fjac[4, 4, i, k] = ( 2.0e+00 - c2 )*
                      u[4, i, j, k, c] * tmp1
               fjac[4, 5, i, k] = c2

               fjac[5, 1, i, k] = ( c2 * 2.0e0 * qs[i, j, k, c]-
                     c1 * ( u[5, i, j, k, c] * tmp1 ) )*
                     ( u[4, i, j, k, c] * tmp1 )
               fjac[5, 2, i, k] = - c2 * ( u[2, i, j, k, c]*u[4, i, j, k, c] )*
                     tmp2
               fjac[5, 3, i, k] = - c2 * ( u[3, i, j, k, c]*u[4, i, j, k, c] )*
                     tmp2
               fjac[5, 4, i, k] = c1 * ( u[5, i, j, k, c] * tmp1 )-
                     c2 * ( qs[i, j, k, c]+
                     u[4, i, j, k, c]*u[4, i, j, k, c] * tmp2 )
               fjac[5, 5, i, k] = c1 * u[4, i, j, k, c] * tmp1

               njac[1, 1, i, k] = 0.0e+00
               njac[1, 2, i, k] = 0.0e+00
               njac[1, 3, i, k] = 0.0e+00
               njac[1, 4, i, k] = 0.0e+00
               njac[1, 5, i, k] = 0.0e+00

               njac[2, 1, i, k] = - c3c4 * tmp2 * u[2, i, j, k, c]
               njac[2, 2, i, k] =   c3c4 * tmp1
               njac[2, 3, i, k] =   0.0e+00
               njac[2, 4, i, k] =   0.0e+00
               njac[2, 5, i, k] =   0.0e+00

               njac[3, 1, i, k] = - c3c4 * tmp2 * u[3, i, j, k, c]
               njac[3, 2, i, k] =   0.0e+00
               njac[3, 3, i, k] =   c3c4 * tmp1
               njac[3, 4, i, k] =   0.0e+00
               njac[3, 5, i, k] =   0.0e+00

               njac[4, 1, i, k] = - con43 * c3c4 * tmp2 * u[4, i, j, k, c]
               njac[4, 2, i, k] =   0.0e+00
               njac[4, 3, i, k] =   0.0e+00
               njac[4, 4, i, k] =   con43 * c3 * c4 * tmp1
               njac[4, 5, i, k] =   0.0e+00

               njac[5, 1, i, k] = - (  c3c4-
                     c1345 ) * tmp3 * (u[2, i, j, k, c]^2)-
                     ( c3c4 - c1345 ) * tmp3 * (u[3, i, j, k, c]^2)-
                     ( con43 * c3c4-
                     c1345 ) * tmp3 * (u[4, i, j, k, c]^2)-
                     c1345 * tmp2 * u[5, i, j, k, c]

               njac[5, 2, i, k] = (  c3c4 - c1345 ) * tmp2 * u[2, i, j, k, c]
               njac[5, 3, i, k] = (  c3c4 - c1345 ) * tmp2 * u[3, i, j, k, c]
               njac[5, 4, i, k] = ( con43 * c3c4-
                     c1345 ) * tmp2 * u[4, i, j, k, c]
               njac[5, 5, i, k] = ( c1345 )* tmp1


            end
         end

#---------------------------------------------------------------------
#     now joacobians set, so form left hand side in z direction
#---------------------------------------------------------------------
         for k = cell_start[3, c]:ksize-cell_end[3, c]
            for i = cell_start[1, c]:isize

               tmp1 = dt * tz1
               tmp2 = dt * tz2

                lhsa[1, 1, i, k] = - tmp2 * fjac[1, 1, i, k-1]-
                     tmp1 * njac[1, 1, i, k-1]-
                     tmp1 * dz1
                lhsa[1, 2, i, k] = - tmp2 * fjac[1, 2, i, k-1]-
                     tmp1 * njac[1, 2, i, k-1]
                lhsa[1, 3, i, k] = - tmp2 * fjac[1, 3, i, k-1]-
                     tmp1 * njac[1, 3, i, k-1]
                lhsa[1, 4, i, k] = - tmp2 * fjac[1, 4, i, k-1]-
                     tmp1 * njac[1, 4, i, k-1]
                lhsa[1, 5, i, k] = - tmp2 * fjac[1, 5, i, k-1]-
                     tmp1 * njac[1, 5, i, k-1]

                lhsa[2, 1, i, k] = - tmp2 * fjac[2, 1, i, k-1]-
                     tmp1 * njac[2, 1, i, k-1]
                lhsa[2, 2, i, k] = - tmp2 * fjac[2, 2, i, k-1]-
                     tmp1 * njac[2, 2, i, k-1]-
                     tmp1 * dz2
                lhsa[2, 3, i, k] = - tmp2 * fjac[2, 3, i, k-1]-
                     tmp1 * njac[2, 3, i, k-1]
                lhsa[2, 4, i, k] = - tmp2 * fjac[2, 4, i, k-1]-
                     tmp1 * njac[2, 4, i, k-1]
                lhsa[2, 5, i, k] = - tmp2 * fjac[2, 5, i, k-1]-
                     tmp1 * njac[2, 5, i, k-1]

                lhsa[3, 1, i, k] = - tmp2 * fjac[3, 1, i, k-1]-
                     tmp1 * njac[3, 1, i, k-1]
                lhsa[3, 2, i, k] = - tmp2 * fjac[3, 2, i, k-1]-
                     tmp1 * njac[3, 2, i, k-1]
                lhsa[3, 3, i, k] = - tmp2 * fjac[3, 3, i, k-1]-
                     tmp1 * njac[3, 3, i, k-1]-
                     tmp1 * dz3
                lhsa[3, 4, i, k] = - tmp2 * fjac[3, 4, i, k-1]-
                     tmp1 * njac[3, 4, i, k-1]
                lhsa[3, 5, i, k] = - tmp2 * fjac[3, 5, i, k-1]-
                     tmp1 * njac[3, 5, i, k-1]

                lhsa[4, 1, i, k] = - tmp2 * fjac[4, 1, i, k-1]-
                     tmp1 * njac[4, 1, i, k-1]
                lhsa[4, 2, i, k] = - tmp2 * fjac[4, 2, i, k-1]-
                     tmp1 * njac[4, 2, i, k-1]
                lhsa[4, 3, i, k] = - tmp2 * fjac[4, 3, i, k-1]-
                     tmp1 * njac[4, 3, i, k-1]
                lhsa[4, 4, i, k] = - tmp2 * fjac[4, 4, i, k-1]-
                     tmp1 * njac[4, 4, i, k-1]-
                     tmp1 * dz4
                lhsa[4, 5, i, k] = - tmp2 * fjac[4, 5, i, k-1]-
                     tmp1 * njac[4, 5, i, k-1]

                lhsa[5, 1, i, k] = - tmp2 * fjac[5, 1, i, k-1]-
                     tmp1 * njac[5, 1, i, k-1]
                lhsa[5, 2, i, k] = - tmp2 * fjac[5, 2, i, k-1]-
                     tmp1 * njac[5, 2, i, k-1]
                lhsa[5, 3, i, k] = - tmp2 * fjac[5, 3, i, k-1]-
                     tmp1 * njac[5, 3, i, k-1]
                lhsa[5, 4, i, k] = - tmp2 * fjac[5, 4, i, k-1]-
                     tmp1 * njac[5, 4, i, k-1]
                lhsa[5, 5, i, k] = - tmp2 * fjac[5, 5, i, k-1]-
                     tmp1 * njac[5, 5, i, k-1]-
                     tmp1 * dz5

                lhsb[1, 1, i, k] = 1.0e+00+
                     tmp1 * 2.0e+00 * njac[1, 1, i, k]+
                     tmp1 * 2.0e+00 * dz1
                lhsb[1, 2, i, k] = tmp1 * 2.0e+00 * njac[1, 2, i, k]
                lhsb[1, 3, i, k] = tmp1 * 2.0e+00 * njac[1, 3, i, k]
                lhsb[1, 4, i, k] = tmp1 * 2.0e+00 * njac[1, 4, i, k]
                lhsb[1, 5, i, k] = tmp1 * 2.0e+00 * njac[1, 5, i, k]

                lhsb[2, 1, i, k] = tmp1 * 2.0e+00 * njac[2, 1, i, k]
                lhsb[2, 2, i, k] = 1.0e+00+
                     tmp1 * 2.0e+00 * njac[2, 2, i, k]+
                     tmp1 * 2.0e+00 * dz2
                lhsb[2, 3, i, k] = tmp1 * 2.0e+00 * njac[2, 3, i, k]
                lhsb[2, 4, i, k] = tmp1 * 2.0e+00 * njac[2, 4, i, k]
                lhsb[2, 5, i, k] = tmp1 * 2.0e+00 * njac[2, 5, i, k]

                lhsb[3, 1, i, k] = tmp1 * 2.0e+00 * njac[3, 1, i, k]
                lhsb[3, 2, i, k] = tmp1 * 2.0e+00 * njac[3, 2, i, k]
                lhsb[3, 3, i, k] = 1.0e+00+
                     tmp1 * 2.0e+00 * njac[3, 3, i, k]+
                     tmp1 * 2.0e+00 * dz3
                lhsb[3, 4, i, k] = tmp1 * 2.0e+00 * njac[3, 4, i, k]
                lhsb[3, 5, i, k] = tmp1 * 2.0e+00 * njac[3, 5, i, k]

                lhsb[4, 1, i, k] = tmp1 * 2.0e+00 * njac[4, 1, i, k]
                lhsb[4, 2, i, k] = tmp1 * 2.0e+00 * njac[4, 2, i, k]
                lhsb[4, 3, i, k] = tmp1 * 2.0e+00 * njac[4, 3, i, k]
                lhsb[4, 4, i, k] = 1.0e+00+
                     tmp1 * 2.0e+00 * njac[4, 4, i, k]+
                     tmp1 * 2.0e+00 * dz4
                lhsb[4, 5, i, k] = tmp1 * 2.0e+00 * njac[4, 5, i, k]

                lhsb[5, 1, i, k] = tmp1 * 2.0e+00 * njac[5, 1, i, k]
                lhsb[5, 2, i, k] = tmp1 * 2.0e+00 * njac[5, 2, i, k]
                lhsb[5, 3, i, k] = tmp1 * 2.0e+00 * njac[5, 3, i, k]
                lhsb[5, 4, i, k] = tmp1 * 2.0e+00 * njac[5, 4, i, k]
                lhsb[5, 5, i, k] = 1.0e+00+
                     tmp1 * 2.0e+00 * njac[5, 5, i, k]+
                     tmp1 * 2.0e+00 * dz5

                lhsc[1, 1, i, j, k, c] =  tmp2 * fjac[1, 1, i, k+1]-
                     tmp1 * njac[1, 1, i, k+1]-
                     tmp1 * dz1
                lhsc[1, 2, i, j, k, c] =  tmp2 * fjac[1, 2, i, k+1]-
                     tmp1 * njac[1, 2, i, k+1]
                lhsc[1, 3, i, j, k, c] =  tmp2 * fjac[1, 3, i, k+1]-
                     tmp1 * njac[1, 3, i, k+1]
                lhsc[1, 4, i, j, k, c] =  tmp2 * fjac[1, 4, i, k+1]-
                     tmp1 * njac[1, 4, i, k+1]
                lhsc[1, 5, i, j, k, c] =  tmp2 * fjac[1, 5, i, k+1]-
                     tmp1 * njac[1, 5, i, k+1]

                lhsc[2, 1, i, j, k, c] =  tmp2 * fjac[2, 1, i, k+1]-
                     tmp1 * njac[2, 1, i, k+1]
                lhsc[2, 2, i, j, k, c] =  tmp2 * fjac[2, 2, i, k+1]-
                     tmp1 * njac[2, 2, i, k+1]-
                     tmp1 * dz2
                lhsc[2, 3, i, j, k, c] =  tmp2 * fjac[2, 3, i, k+1]-
                     tmp1 * njac[2, 3, i, k+1]
                lhsc[2, 4, i, j, k, c] =  tmp2 * fjac[2, 4, i, k+1]-
                     tmp1 * njac[2, 4, i, k+1]
                lhsc[2, 5, i, j, k, c] =  tmp2 * fjac[2, 5, i, k+1]-
                     tmp1 * njac[2, 5, i, k+1]

                lhsc[3, 1, i, j, k, c] =  tmp2 * fjac[3, 1, i, k+1]-
                     tmp1 * njac[3, 1, i, k+1]
                lhsc[3, 2, i, j, k, c] =  tmp2 * fjac[3, 2, i, k+1]-
                     tmp1 * njac[3, 2, i, k+1]
                lhsc[3, 3, i, j, k, c] =  tmp2 * fjac[3, 3, i, k+1]-
                     tmp1 * njac[3, 3, i, k+1]-
                     tmp1 * dz3
                lhsc[3, 4, i, j, k, c] =  tmp2 * fjac[3, 4, i, k+1]-
                     tmp1 * njac[3, 4, i, k+1]
                lhsc[3, 5, i, j, k, c] =  tmp2 * fjac[3, 5, i, k+1]-
                     tmp1 * njac[3, 5, i, k+1]

                lhsc[4, 1, i, j, k, c] =  tmp2 * fjac[4, 1, i, k+1]-
                     tmp1 * njac[4, 1, i, k+1]
                lhsc[4, 2, i, j, k, c] =  tmp2 * fjac[4, 2, i, k+1]-
                     tmp1 * njac[4, 2, i, k+1]
                lhsc[4, 3, i, j, k, c] =  tmp2 * fjac[4, 3, i, k+1]-
                     tmp1 * njac[4, 3, i, k+1]
                lhsc[4, 4, i, j, k, c] =  tmp2 * fjac[4, 4, i, k+1]-
                     tmp1 * njac[4, 4, i, k+1]-
                     tmp1 * dz4
                lhsc[4, 5, i, j, k, c] =  tmp2 * fjac[4, 5, i, k+1]-
                     tmp1 * njac[4, 5, i, k+1]

                lhsc[5, 1, i, j, k, c] =  tmp2 * fjac[5, 1, i, k+1]-
                     tmp1 * njac[5, 1, i, k+1]
                lhsc[5, 2, i, j, k, c] =  tmp2 * fjac[5, 2, i, k+1]-
                     tmp1 * njac[5, 2, i, k+1]
                lhsc[5, 3, i, j, k, c] =  tmp2 * fjac[5, 3, i, k+1]-
                     tmp1 * njac[5, 3, i, k+1]
                lhsc[5, 4, i, j, k, c] =  tmp2 * fjac[5, 4, i, k+1]-
                     tmp1 * njac[5, 4, i, k+1]
                lhsc[5, 5, i, j, k, c] =  tmp2 * fjac[5, 5, i, k+1]-
                     tmp1 * njac[5, 5, i, k+1]-
                     tmp1 * dz5

            end
         end


#---------------------------------------------------------------------
#     outer most do loops - sweeping in i direction
#---------------------------------------------------------------------
         if FIRST == 1

#---------------------------------------------------------------------
#     multiply c[i,j,kstart] by b_inverse and copy back to !
#     multiply rhs[kstart] by b_inverse(kstart) and copy to rhs
#---------------------------------------------------------------------
#dir$ ivdep
            for i = cell_start[1, c]:isize
               binvcrhs(  lhsb[1, 1, i, kstart],
                               lhsc[1, 1, i, j, kstart, c],
                              rhs[1, i, j, kstart, c] )
            end

         end

#---------------------------------------------------------------------
#     begin inner most do loop
#     do all the elements of the cell unless last 
#---------------------------------------------------------------------
         for k = kstart+FIRST:ksize-LAST
#dir$ ivdep
            for i = cell_start[1, c]:isize

#---------------------------------------------------------------------
#     subtract A*lhs_vector(k-1) from lhs_vector(k)
#     
#     rhs[k] = rhs[k] - A*rhs[k-1]
#---------------------------------------------------------------------
               matvec_sub( lhsa[1, 1, i, k],
                               rhs[1, i, j, k-1, c], rhs[1, i, j, k, c])

#---------------------------------------------------------------------
#      b[k] =  b[k] - C(k-1)*A(k)
#     call matmul_sub(aa,i,j,k,c,cc,i,j,k-1,c,bb,i,j,k,c)
#---------------------------------------------------------------------
               matmul_sub( lhsa[1, 1, i, k],
                                lhsc[1, 1, i, j, k-1, c],
                                lhsb[1, 1, i, k])

#---------------------------------------------------------------------
#     multiply c[i,j,k] by b_inverse and copy back to !
#     multiply rhs[i,j,1] by b_inverse(i,j,1) and copy to rhs
#---------------------------------------------------------------------
               binvcrhs(  lhsb[1, 1, i, k],
                               lhsc[1, 1, i, j, k, c],
                              rhs[1, i, j, k, c] )

            end
         end

#---------------------------------------------------------------------
#     Now finish up special cases for last cell
#---------------------------------------------------------------------
         if LAST == 1

#dir$ ivdep
            for i = cell_start[1, c]:isize
#---------------------------------------------------------------------
#     rhs[ksize] = rhs[ksize] - A*rhs[ksize-1]
#---------------------------------------------------------------------
               matvec_sub( lhsa[1, 1, i, ksize],
                               rhs[1, i, j, ksize-1, c], rhs[1, i, j, ksize, c])

#---------------------------------------------------------------------
#      b[ksize] =  b[ksize] - C(ksize-1)*A(ksize)
#     call matmul_sub(aa,i,j,ksize,c,
#     $              cc,i,j,ksize-1,c,bb,i,j,ksize,c)
#---------------------------------------------------------------------
               matmul_sub( lhsa[1, 1, i, ksize],
                                lhsc[1, 1, i, j, ksize-1, c],
                                lhsb[1, 1, i, ksize])

#---------------------------------------------------------------------
#     multiply rhs[ksize] by b_inverse(ksize) and copy to rhs
#---------------------------------------------------------------------
               binvrhs(  lhsb[1, 1, i, ksize],
                             rhs[1, i, j, ksize, c] )
            end

         end
      end


      return nothing
      end






