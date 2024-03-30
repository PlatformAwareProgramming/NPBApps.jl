#---------------------------------------------------------------------
#
#   set up the sub-domain sizes
#
#---------------------------------------------------------------------

 function subdomain(iz, row, col, west, east, north, south, nx0, ny0, nz0, nx, ny, nz, ipt, jpt, ist, jst, iend, jend)

#---------------------------------------------------------------------
#   x dimension
#---------------------------------------------------------------------
      mm   = mod(nx0, xdim)
      if row <= mm
        nx[iz] = div(nx0,xdim) + 1
        ipt[iz] = (row-1)*nx[iz]
      else
        nx[iz] = div(nx0,xdim)
        ipt[iz] = (row-1)*nx[iz] + mm
      end

#---------------------------------------------------------------------
#   y dimension
#---------------------------------------------------------------------
      mm   = mod(ny0, ydim)
      if col <= mm
        ny[iz] = div(ny0,ydim) + 1
        jpt[iz] = (col-1)*ny[iz]
      else
        ny[iz] = div(ny0,ydim)
        jpt[iz] = (col-1)*ny[iz] + mm
      end

#---------------------------------------------------------------------
#   z dimension
#---------------------------------------------------------------------
      nz[iz] = nz0

#---------------------------------------------------------------------
#   check the sub-domain size
#---------------------------------------------------------------------
      if (nx[iz] < 3) || (ny[iz] < 3) || (nz[iz] < 3)
          @printf(stdout, "     SUBDOMAIN SIZE IS TOO SMALL - \n     ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS\n     SO THAT NX, NY AND NZ ARE GREATER THAN OR EQUAL\n     TO 3 THEY ARE CURRENTLY%5i%5i%5i\n", nx, ny, nz)
          ERRORCODE = 1
          MPI.Abort( MPI.COMM_WORLD, ERRORCODE)
      end

      if (nx[iz] > nx0) ||(ny[iz] > ny0) || (nz[iz] > nz0)
          @printf(stdout, "     SUBDOMAIN SIZE IS TOO LARGE - \n     ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS\n     SO THAT NX, NY AND NZ ARE LESS THAN OR EQUAL TO \n     ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY.  THEY ARE\n     CURRENTLY%5i%5i%5i\n", nx, ny, nz)
          ERRORCODE = 1
          MPI.Abort( MPI.COMM_WORLD, ERRORCODE)
      end


#---------------------------------------------------------------------
#   set up the start and end in i and j extents for all processors
#---------------------------------------------------------------------
      ist[iz] = 1
      iend[iz] = nx[iz]
      if (north == -1) ist[iz] = 2 end
      if (south == -1) iend[iz] = nx[iz] - 1 end

      jst[iz] = 1
      jend[iz] = ny[iz]
      if (west == -1) jst[iz] = 2 end
      if (east == -1) jend[iz] = ny[iz] - 1 end

      return nothing
end


