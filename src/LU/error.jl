#---------------------------------------------------------------------
#
#   compute the solution error
#
#---------------------------------------------------------------------

 function ERROR(u, nx, ny, nz, ipt, jpt, ist, jst, iend, jend, errnm, comm_solve)

      u000ijk = Array{FloatType}(undef, 5)
      dummy = zeros(FloatType, 5)

      errnm .= 0.0e+00

      for k = 2:nz-1
         for j = jst:jend
            jglob = jpt + j
            for i = ist:iend
               iglob = ipt + i
               exact( iglob, jglob, k, u000ijk, nx, ny, nz)
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
      # MPI_ALLREDUCE( dummy, errnm, 5, dp_type, MPI_SUM, comm_solve, IERROR )
       MPI.Allreduce!(dummy, errnm, MPI.SUM, comm_solve)

      for m = 1:5
         errnm[m] = sqrt( errnm[m] / ( float(nx-2)*(ny-2)*(nz-2) ) )
      end

#      if (id.eq.0) then
#        write (*,1002) ( errnm(m), m = 1, 5 )
#      end if

# 1002 format(1x/1x,'RMS-norm of error in soln. to ',  		       'first pde  = ',1pe12.5/,  		       1x,'RMS-norm of error in soln. to ',  		       'second pde = ',1pe12.5/,  		       1x,'RMS-norm of error in soln. to ',  		       'third pde  = ',1pe12.5/,  		       1x,'RMS-norm of error in soln. to ',  		       'fourth pde = ',1pe12.5/,  		       1x,'RMS-norm of error in soln. to ',  		       'fifth pde  = ',1pe12.5)

      return nothing
end
