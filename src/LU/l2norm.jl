using FortranFiles
using OffsetArrays
using Parameters
using Printf


#---------------------------------------------------------------------
#---------------------------------------------------------------------
      function l2norm( ldx, ldy, ldz,
                          nx0, ny0, nz0,
                          ist, iend,
                          jst, jend,
                          v, SUM )
#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#   to compute the l2-norm of vector v.
#---------------------------------------------------------------------

#      use timing
#      use mpinpb

#      implicit none

#---------------------------------------------------------------------
#  input parameters
#---------------------------------------------------------------------
#      integer ldx, ldy, ldz
#      integer nx0, ny0, nz0
#      integer ist, iend
#      integer jst, jend
#      DOUBLEPRECISION  v[5,-1:ldx+2,-1:ldy+2,*], SUM[5]

#---------------------------------------------------------------------
#  local variables
#---------------------------------------------------------------------
#      integer i, j, k, m
#      DOUBLEPRECISION  dummy[5]

#      integer IERROR


      for m = 1:5
         dummy[m] = 0.0e+00
      end

      for k = 2:nz0-1
         for j = jst:jend
            for i = ist:iend
               for m = 1:5
                  dummy[m] = dummy[m] + v[m,i,j,k] * v[m,i,j,k]
               end
            end
         end
      end

#---------------------------------------------------------------------
#   compute the global sum of individual contributions to dot product.
#---------------------------------------------------------------------
      if (timeron) timer_start(t_rcomm) end
      MPI_ALLREDUCE( dummy,
                          SUM,
                          5,
                          dp_type,
                          MPI_SUM,
                          comm_solve,
                          IERROR )
      if (timeron) timer_stop(t_rcomm) end

      for m = 1:5
         SUM[m] = sqrt( SUM[m] / ( float(nx0-2)*(ny0-2)*(nz0-2) ) )
      end

      return nothing
      end
