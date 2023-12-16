using FortranFiles
using OffsetArrays
using Parameters
using Printf


#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function subdomain()

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#      use lu_data
#      use mpinpb

#      implicit none

#---------------------------------------------------------------------
#  local variables
#---------------------------------------------------------------------
#      integer mm, ierror, errorcode


#---------------------------------------------------------------------
#
#   set up the sub-domain sizes
#
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#   x dimension
#---------------------------------------------------------------------
      mm   = mod(nx0, xdim)
      if row <= mm
        nx = nx0/xdim + 1
        ipt = (row-1)*nx
      else
        nx = nx0/xdim
        ipt = (row-1)*nx + mm
      end

#---------------------------------------------------------------------
#   y dimension
#---------------------------------------------------------------------
      mm   = mod(ny0, ydim)
      if col <= mm
        ny = ny0/ydim + 1
        jpt = (col-1)*ny
      else
        ny = ny0/ydim
        jpt = (col-1)*ny + mm
      end

#---------------------------------------------------------------------
#   z dimension
#---------------------------------------------------------------------
      nz = nz0

#---------------------------------------------------------------------
#   check the sub-domain size
#---------------------------------------------------------------------
      if ( nx < 3 ) ||(
            ny < 3 ) ||(
            nz < 3 )
         @printf(stdout, "     SUBDOMAIN SIZE IS TOO SMALL - \n     ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS\n     SO THAT NX, NY AND NZ ARE GREATER THAN OR EQUAL\n     TO 3 THEY ARE CURRENTLY%5i%5i%5i\n", nx, ny, nz)
# 2001    format(5x,'SUBDOMAIN SIZE IS TOO SMALL - ',  		              /5x,'ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS',  		              /5x,'SO THAT NX, NY AND NZ ARE GREATER THAN OR EQUAL',  		              /5x,'TO 3 THEY ARE CURRENTLY', 3I5)
          ERRORCODE = 1
          MPI_ABORT( MPI_COMM_WORLD,
                          ERRORCODE,
                          IERROR )
      end

      if ( nx > isiz1 ) ||(
            ny > isiz2 ) ||(
            nz > isiz3 )
         @printf(stdout, "     SUBDOMAIN SIZE IS TOO LARGE - \n     ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS\n     SO THAT NX, NY AND NZ ARE LESS THAN OR EQUAL TO \n     ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY.  THEY ARE\n     CURRENTLY%5i%5i%5i\n", nx, ny, nz)
# 2002    format(5x,'SUBDOMAIN SIZE IS TOO LARGE - ',  		              /5x,'ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS',  		              /5x,'SO THAT NX, NY AND NZ ARE LESS THAN OR EQUAL TO ',  		              /5x,'ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY.  THEY ARE',  		              /5x,'CURRENTLY', 3I5)
          ERRORCODE = 1
          MPI_ABORT( MPI_COMM_WORLD,
                          ERRORCODE,
                          IERROR )
      end


#---------------------------------------------------------------------
#   set up the start and end in i and j extents for all processors
#---------------------------------------------------------------------
      ist = 1
      iend = nx
      if (north == -1) ist = 2 end
      if (south == -1) iend = nx - 1 end

      jst = 1
      jend = ny
      if (west == -1) jst = 2 end
      if (east == -1) jend = ny - 1 end

      return nothing
      end


