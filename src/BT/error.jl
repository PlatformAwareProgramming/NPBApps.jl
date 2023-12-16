using FortranFiles
using OffsetArrays
using Parameters
using Printf

#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function error_norm(rms)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     this function computes the norm of the difference between the
#     computed solution and the exact solution
#---------------------------------------------------------------------

#      use bt_data
#      use mpinpb

#      implicit none

#      integer c, i, j, k, m, ii, jj, kk, d, ERROR
#      DOUBLEPRECISION xi, eta, zeta, u_exact[5], rms[5], rms_work[5],  
#           add

      timer_start(t_enorm)

      for m = 1:5
         rms_work[m] = 0.0e0
      end

      for c = 1:ncells
         kk = 0
         for k = cell_low[3, c]:cell_high[3, c]
            zeta = float(k) * dnzm1
            jj = 0
            for j = cell_low[2, c]:cell_high[2, c]
               eta = float(j) * dnym1
               ii = 0
               for i = cell_low[1, c]:cell_high[1, c]
                  xi = float(i) * dnxm1
                  exact_solution(xi, eta, zeta, u_exact)

                  for m = 1:5
                     add = u[m, ii, jj, kk, c]-u_exact[m]
                     rms_work[m] = rms_work[m] + add*add
                  end
                  ii = ii + 1
               end
               jj = jj + 1
            end
            kk = kk + 1
         end
      end

      mpi_allreduce(rms_work, rms, 5, dp_type,
           MPI_SUM, comm_setup, ERROR)

      for m = 1:5
         for d = 1:3
            rms[m] = rms[m] / float(predecessor[d]-2)
         end
         rms[m] = sqrt(rms[m])
      end

      timer_stop(t_enorm)

      return nothing
      end


#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function rhs_norm(rms)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#      use bt_data
#      use mpinpb

#      implicit none

#      integer c, i, j, k, d, m, ERROR
#      DOUBLEPRECISION rms[5], rms_work[5], add

      for m = 1:5
         rms_work[m] = 0.0e0
      end

      for c = 1:ncells
         for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
            for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
               for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                  for m = 1:5
                     add = rhs[m, i, j, k, c]
                     rms_work[m] = rms_work[m] + add*add
                  end
               end
            end
         end
      end

      mpi_allreduce(rms_work, rms, 5, dp_type,
           MPI_SUM, comm_setup, ERROR)

      for m = 1:5
         for d = 1:3
            rms[m] = rms[m] / float(predecessor[d]-2)
         end
         rms[m] = sqrt(rms[m])
      end

      return nothing
      end

