using FortranFiles
using OffsetArrays
using Parameters
using Printf


#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function isqrt2(i, xdim)

#      implicit none
#      integer i, xdim

#      integer ydim, square

      xdim = -1
      if (i <= 0) return end

      square = 0;
      xdim = 1
      while (square <= i)
         square = xdim*xdim
         if (square == i) return end
         xdim = xdim + 1
      end

      xdim = xdim - 1
      ydim = i / xdim
      while (xdim*ydim != i && 2*ydim >= xdim)
         xdim = xdim + 1
         ydim = i / xdim
      end

      if (xdim*ydim != i || 2*ydim < xdim) xdim = -1 end

      return nothing
      end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function proc_grid()

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#      use lu_data
#      use mpinpb

#      implicit none

#---------------------------------------------------------------------
#  local variables
#---------------------------------------------------------------------
#      integer xdim0, ydim0, IERROR
#      integer xdiv, ydiv

#---------------------------------------------------------------------
#  calculate sub-domain array size
#---------------------------------------------------------------------
      isqrt2(no_nodes, xdiv)

      if xdiv <= 0
         if (id == 0) @printf(stdout, " ERROR: could not determine proper proc_grid for nprocs = %6i\n", no_nodes) end
         @label L2000
         format(" ERROR: could not determine proper proc_grid for nprocs = ", i6)
         MPI_ABORT( MPI_COMM_WORLD, MPI_ERR_OTHER, IERROR )
         exit(1)
      end

      ydiv = no_nodes/xdiv
      isiz1 = isiz01/xdiv
      if (isiz1*xdiv < isiz01) isiz1 = isiz1 + 1 end
      isiz2 = isiz02/ydiv
      if (isiz2*ydiv < isiz02) isiz2 = isiz2 + 1 end
      nnodes_xdim = xdiv
      isiz3 = isiz03

#---------------------------------------------------------------------
#
#   set up a two-d grid for processors: column-major ordering of unknowns
#
#---------------------------------------------------------------------

      xdim0  = nnodes_xdim
      ydim0  = no_nodes/xdim0

      ydim   = sqrt(float(num))+0.001e0
      xdim   = num/ydim
      while (ydim >= ydim0 && xdim*ydim != num)
         ydim = ydim - 1
         xdim = num/ydim
      end

      if xdim < xdim0 || ydim < ydim0 ||
          xdim*ydim != num
         if (id == 0) @printf(stdout, " ERROR: could not determine proper proc_grid for nprocs = %6i\n", num) end
         MPI_ABORT( MPI_COMM_WORLD, MPI_ERR_OTHER, IERROR )
         exit(1)
      end

      if (id == 0 && num != 2^ndim)
         @printf(stdout, " Proc_grid for nprocs =%6i:%5i x%5i\n\n", num, xdim, ydim)
      end
      @label L2100
      format(" Proc_grid for nprocs =", i6, ':', i5, " x",i5)

      row    = mod(id, xdim) + 1
      col    = id/xdim + 1


      return nothing
      end


