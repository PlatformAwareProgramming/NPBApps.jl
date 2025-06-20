#---------------------------------------------------------------------
#   to compute the l2-norm of vector v.
#---------------------------------------------------------------------

dummy = zeros(Float64, 5)

function l2norm( ldx, ldy, ldz, nx0, ny0, nz0, ist, iend, jst, jend, v, SUM, 
   comm_solve, 
   timeron,
   )

      dummy .= 0

      @inbounds for k = 2:nz0-1
         for j = jst:jend
            for i = ist:iend
               for m = 1:5
                  dummy[m] += v[m,i,j,k] * v[m,i,j,k]
               end
            end
         end
      end

#---------------------------------------------------------------------
#   compute the global sum of individual contributions to dot product.
#---------------------------------------------------------------------
     if timeron timer_start(t_rcomm) end

      MPI.Allreduce!(dummy, SUM, MPI.SUM, comm_solve)

#      MPI_ALLREDUCE( dummy, SUM, 5, dp_type, MPI_SUM, comm_solve, IERROR )
     if timeron timer_stop(t_rcomm) end

     @inbounds for m = 1:5
         SUM[m] = sqrt(SUM[m] / (float(nx0-2)*(ny0-2)*(nz0-2)))
      end 

      return nothing
end
