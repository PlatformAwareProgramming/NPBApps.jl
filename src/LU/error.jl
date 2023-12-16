using FortranFiles
using OffsetArrays
using Parameters
using Printf

#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function ERROR()

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#
#   compute the solution error
#
#---------------------------------------------------------------------

#      use lu_data
#      use mpinpb

#      implicit none

#---------------------------------------------------------------------
#  local variables
#---------------------------------------------------------------------
#      integer i, j, k, m
#      integer iglob, jglob
#      DOUBLEPRECISION  tmp
#      DOUBLEPRECISION  u000ijk[5], dummy[5]

#      integer IERROR


      for m = 1:5
         errnm(m) = 0.0e+00
         dummy[m] = 0.0e+00
      end

      for k = 2:nz-1
         for j = jst:jend
            jglob = jpt + j
            for i = ist:iend
               iglob = ipt + i
               exact( iglob, jglob, k, u000ijk )
               for m = 1:5
                  tmp = ( u000ijk[m] - u[m, i, j, k] )
                  dummy[m] = dummy[m] + tmp ^ 2
               end
            end
         end
      end

#---------------------------------------------------------------------
#   compute the global sum of individual contributions to dot product.
#---------------------------------------------------------------------
      MPI_ALLREDUCE( dummy,
                          errnm,
                          5,
                          dp_type,
                          MPI_SUM,
                          comm_solve,
                          IERROR )

      for m = 1:5
         errnm(m) = sqrt( errnm(m) / ( float(nx0-2)*(ny0-2)*(nz0-2) ) )
      end

#      if (id.eq.0) then
#        write (*,1002) ( errnm(m), m = 1, 5 )
#      end if

# 1002 format(1x/1x,'RMS-norm of error in soln. to ',  		       'first pde  = ',1pe12.5/,  		       1x,'RMS-norm of error in soln. to ',  		       'second pde = ',1pe12.5/,  		       1x,'RMS-norm of error in soln. to ',  		       'third pde  = ',1pe12.5/,  		       1x,'RMS-norm of error in soln. to ',  		       'fourth pde = ',1pe12.5/,  		       1x,'RMS-norm of error in soln. to ',  		       'fifth pde  = ',1pe12.5)

      return nothing
      end
