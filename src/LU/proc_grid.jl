


#---------------------------------------------------------------------
#---------------------------------------------------------------------

 function isqrt2(i)

      xdim = -1
      if (i <= 0) 
         return xdim 
      end

      square = 0;
      xdim = 1
      while (square <= i)
         square = xdim*xdim
         if (square == i) 
            return xdim 
         end
         xdim = xdim + 1
      end

      xdim = xdim - 1
      ydim = i / xdim
      while (xdim*ydim != i && 2*ydim >= xdim)
         xdim = xdim + 1
         ydim = i / xdim
      end

      if (xdim*ydim != i || 2*ydim < xdim) 
         xdim = -1 
      end

      return xdim
end


#---------------------------------------------------------------------
#  calculate sub-domain array size
#---------------------------------------------------------------------

function proc_grid(isiz01, isiz02, isiz03)

      xdiv = isqrt2(no_nodes)

      if xdiv <= 0
         if (id == 0) @printf(stdout, " ERROR: could not determine proper proc_grid for nprocs = %6i\n", no_nodes) end
         @label L2000
         format(" ERROR: could not determine proper proc_grid for nprocs = ", i6)
         MPI.Abort( MPI.COMM_WORLD, MPI.ERR_OTHER )
         exit(1)
      end

      ydiv = div(no_nodes,xdiv)
      global isiz1 = div(isiz01,xdiv)
      if (isiz1*xdiv < isiz01) isiz1 = isiz1 + 1 end
      global isiz2 = div(isiz02,ydiv)
      if (isiz2*ydiv < isiz02) isiz2 = isiz2 + 1 end
      nnodes_xdim = xdiv
      global isiz3 = isiz03

#---------------------------------------------------------------------
#
#   set up a two-d grid for processors: column-major ordering of unknowns
#
#---------------------------------------------------------------------

      xdim0  = nnodes_xdim
      ydim0  = div(no_nodes,xdim0)

      ydim   = Int(sqrt(num))
      xdim   = div(num,ydim)
      
      while (ydim >= ydim0 && xdim*ydim != num)
         ydim = ydim - 1
         xdim = div(num,ydim)
      end

      if xdim < xdim0 || ydim < ydim0 || xdim*ydim != num
         if (id == 0) @printf(stdout, " ERROR: could not determine proper proc_grid for nprocs = %6i\n", num) end
         MPI.Abort( MPI.COMM_WORLD, MPI.MPI_ERR_OTHER )
         exit(1)
      end

      if (id == 0 && num != 2^ndim)
         @printf(stdout, " Proc_grid for nprocs =%6i:%5i x%5i\n\n", num, xdim, ydim)
      end
      @label L2100

      global row = mod(id,xdim) + 1
      global col = div(id,xdim) + 1

      return nothing
end


