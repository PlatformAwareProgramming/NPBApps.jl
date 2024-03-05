#---------------------------------------------------------------------
# this function computes the norm of the difference between the
# computed solution and the exact solution
#---------------------------------------------------------------------


function error_norm(z, rms, grid_points)

   rms_work = zeros(Float64,5)
   u_exact = Array{Float64}(undef,5)   
   
   for c = 1:ncells
      kk = 0
      for k = cell_low[z][3, c]:cell_high[z][3, c]
         zeta = float(k) * dnzm1
         jj = 0
         for j = cell_low[z][2, c]:cell_high[z][2, c]
            eta = float(j) * dnym1
            ii = 0
            for i = cell_low[z][1, c]:cell_high[z][1, c]
               xi = float(i) * dnxm1
               u_exact = exact_solution(xi, eta, zeta)

               for m = 1:5
                  add = u[z][m, ii, jj, kk, c] - u_exact[m]
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
      rms[m] = rms[m] / (float(grid_points[1]-2)*float(grid_points[2]-2)*float(grid_points[3]-2))
      rms[m] = sqrt(rms[m])
   end

   return nothing
end



function rhs_norm(z, rms, grid_points)

   rms_work = zeros(Float64, 5)

   for c = 1:ncells
      kk = cell_low[z][3, c] + cell_start[z][3, c]
      for k = cell_start[z][3, c]:cell_size[z][3, c]-cell_end[z][3, c]-1
         jj = cell_low[z][2, c] + cell_start[z][2, c]
         for j = cell_start[z][2, c]:cell_size[z][2, c]-cell_end[z][2, c]-1
            ii = cell_low[z][1, c] + cell_start[z][1, c]
            for i = cell_start[z][1, c]:cell_size[z][1, c]-cell_end[z][1, c]-1
               for m = 1:5
                  add = rhs[z][m, i, j, k, c]
                  rms_work[m] = rms_work[m] + add*add
               end
               ii += 1
            end
            jj += 1
         end
         kk += 1
      end
   end

   MPI.Allreduce!(rms_work, rms, MPI.SUM, comm_setup)

   for m = 1:5
      rms[m] = rms[m] / (float(grid_points[1]-2) * float(grid_points[2]-2) * float(grid_points[3]-2))
      rms[m] = sqrt(rms[m])
   end

   return nothing
end


function write_rhs(z)
   zone = proc_zone_id[z]
   for c = 1:ncells
      kk = cell_low[z][3, c] + cell_start[z][3, c]
      for k = cell_start[z][3, c]:cell_size[z][3, c]-cell_end[z][3, c]-1
         jj = cell_low[z][2, c] + cell_start[z][2, c]
         for j = cell_start[z][2, c]:cell_size[z][2, c]-cell_end[z][2, c]-1
            ii = cell_low[z][1, c] + cell_start[z][1, c]
            for i = cell_start[z][1, c]:cell_size[z][1, c]-cell_end[z][1, c]-1
               for m = 1:5
                  @info "rhs[$m, $ii, $jj, $kk][$zone] = $(rhs[z][m, i, j, k, c])"
               end
               ii += 1
            end
            jj += 1
         end
         kk += 1
      end
   end

   return nothing
end

function write_u(z)
   zone = proc_zone_id[z]
   for c = 1:ncells
      kk = cell_low[z][3, c] + cell_start[z][3, c]
      for k = cell_start[z][3, c]:cell_size[z][3, c]-cell_end[z][3, c]-1
         jj = cell_low[z][2, c] + cell_start[z][2, c]
         for j = cell_start[z][2, c]:cell_size[z][2, c]-cell_end[z][2, c]-1
            ii = cell_low[z][1, c] + cell_start[z][1, c]
            for i = cell_start[z][1, c]:cell_size[z][1, c]-cell_end[z][1, c]-1
               for m = 1:5
                  @info "$node: u[$m, $ii, $jj, $kk][$zone] = $(u[z][m, i, j, k, c])"
               end
               ii += 1
            end
            jj += 1
         end
         kk += 1
      end
   end

   return nothing
end