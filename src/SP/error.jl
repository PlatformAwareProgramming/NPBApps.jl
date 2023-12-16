
rms_work = Array{Float64}(undef,5)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

function error_norm(rms)

#---------------------------------------------------------------------
#---------------------------------------------------------------------



#---------------------------------------------------------------------
# this function computes the norm of the difference between the
# computed solution and the exact solution
#---------------------------------------------------------------------

#       use sp_data
#       use mpinpb

#       implicit none

#       integer c, i, j, k, m, ii, jj, kk, d, ERROR
#       DOUBLEPRECISION xi, eta, zeta, u_exact[5], rms[5], rms_work[5],  
#                        add

      u_exact = Array{Float64}(undef,5)

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
                   u_exact = exact_solution(xi, eta, zeta)

                   for m = 1:5
                      add = u[ii, jj, kk, m, c]-u_exact[m]
                      rms_work[m] = rms_work[m] + add*add
                   end
                   ii = ii + 1
                end
                jj = jj + 1
             end
             kk = kk + 1
          end
       end

#      mpi_allreduce(rms_work, rms, 5, dp_type, MPI_SUM, comm_setup, ERROR)
       MPI.Allreduce!(rms_work, rms, MPI.SUM, comm_setup)

       for m = 1:5
          for d = 1:3
             rms[m] = rms[m] / float(grid_points[d]-2)
          end
          rms[m] = sqrt(rms[m])
       end

       return nothing
       end



       function rhs_norm(rms)

#       use sp_data
#       use mpinpb

#       implicit none

#       integer c, i, j, k, d, m, ERROR
#       DOUBLEPRECISION rms[5], rms_work[5], add

       for m = 1:5
          rms_work[m] = 0.0e0
       end

       for c = 1:ncells
          for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
             for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
                for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                   for m = 1:5
                      add = rhs[i, j, k, m, c]
                      rms_work[m] = rms_work[m] + add*add
                   end
                end
             end
          end
       end



       #mpi_allreduce(rms_work, rms, 5, dp_type, MPI_SUM, comm_setup, ERROR)
       MPI.Allreduce!(rms_work, rms, MPI.SUM, comm_setup)

       for m = 1:5
          for d = 1:3
             rms[m] = rms[m] / float(grid_points[d]-2)
          end
          rms[m] = sqrt(rms[m])
       end

       return nothing
       end


