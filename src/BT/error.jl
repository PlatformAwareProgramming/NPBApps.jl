
#---------------------------------------------------------------------
#     this function computes the norm of the difference between the
#     computed solution and the exact solution
#---------------------------------------------------------------------

function error_norm(ncells, u, rms, grid_points, cell_low, cell_high, dnxm1, dnym1, dnzm1, comm_setup)

      rms_work = zeros(Float64,5)
      u_exact = Array{Float64}(undef,5)

      timer_start(t_enorm)

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

      MPI.Allreduce!(rms_work, rms, MPI.SUM, comm_setup)

      for m = 1:5
         for d = 1:3
            rms[m] = rms[m] / float(grid_points[d]-2)
         end
         rms[m] = sqrt(rms[m])
      end

      timer_stop(t_enorm)

      return nothing
end


#---------------------------------------------------------------------
#---------------------------------------------------------------------

 function rhs_norm(ncells, rhs, rms, grid_points, cell_start, cell_end, cell_size, comm_setup)

      rms_work = zeros(Float64, 5)

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

      MPI.Allreduce!(rms_work, rms, MPI.SUM, comm_setup)

      for m = 1:5
         for d = 1:3
            rms[m] = rms[m] / float(grid_points[d]-2)
         end
         rms[m] = sqrt(rms[m])
      end

      return nothing
end

