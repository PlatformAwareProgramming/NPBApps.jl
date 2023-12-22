#---------------------------------------------------------------------
#
#   set up the sub-domain sizes
#
#---------------------------------------------------------------------

 function subdomain()

#---------------------------------------------------------------------
#   x dimension
#---------------------------------------------------------------------
      mm   = mod(nx0, xdim)
      if row <= mm
        global nx = div(nx0,xdim) + 1
        global ipt = (row-1)*nx
      else
        global nx = div(nx0,xdim)
        global ipt = (row-1)*nx + mm
      end

#---------------------------------------------------------------------
#   y dimension
#---------------------------------------------------------------------
      mm   = mod(ny0, ydim)
      if col <= mm
        global ny = div(ny0,ydim) + 1
        global jpt = (col-1)*ny
      else
        global ny = div(ny0,ydim)
        global jpt = (col-1)*ny + mm
      end

#---------------------------------------------------------------------
#   z dimension
#---------------------------------------------------------------------
      global nz = nz0

#---------------------------------------------------------------------
#   check the sub-domain size
#---------------------------------------------------------------------
      if (nx < 3) || (ny < 3) || (nz < 3)
          @printf(stdout, "     SUBDOMAIN SIZE IS TOO SMALL - \n     ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS\n     SO THAT NX, NY AND NZ ARE GREATER THAN OR EQUAL\n     TO 3 THEY ARE CURRENTLY%5i%5i%5i\n", nx, ny, nz)
          ERRORCODE = 1
          MPI.Abort( MPI.COMM_WORLD, ERRORCODE)
      end

      if (nx > isiz1) ||(ny > isiz2) || (nz > isiz3)
          @printf(stdout, "     SUBDOMAIN SIZE IS TOO LARGE - \n     ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS\n     SO THAT NX, NY AND NZ ARE LESS THAN OR EQUAL TO \n     ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY.  THEY ARE\n     CURRENTLY%5i%5i%5i\n", nx, ny, nz)
          ERRORCODE = 1
          MPI.Abort( MPI.COMM_WORLD, ERRORCODE)
      end


#---------------------------------------------------------------------
#   set up the start and end in i and j extents for all processors
#---------------------------------------------------------------------
      global ist = 1
      global iend = nx
      if (north == -1) ist = 2 end
      if (south == -1) iend = nx - 1 end

      global jst = 1
      global jend = ny
      if (west == -1) jst = 2 end
      if (east == -1) jend = ny - 1 end

      return nothing
end


