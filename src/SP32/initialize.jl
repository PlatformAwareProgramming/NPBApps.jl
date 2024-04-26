#---------------------------------------------------------------------
# This subroutine initializes the field variable u using 
# tri-linear transfinite interpolation of the boundary values     
#---------------------------------------------------------------------

function initialize(z)

      Pface1 = Array{Array{Float32}}(undef,2)
      Pface2 = Array{Array{Float32}}(undef,2)
      Pface3 = Array{Array{Float32}}(undef,2)
      temp = Array{Float32}(undef, 5)

#---------------------------------------------------------------------
#  Later (in compute_rhs) we compute 1/u for every element. A few of 
#  the corner elements are not used, but it convenient (and faster) 
#  to compute the whole thing with a simple loop. Make sure those 
#  values are nonzero by initializing the whole thing here. 
#---------------------------------------------------------------------
      for c = 1:ncells
         for kk = -1:KMAX[z]
            for jj = -1:JMAX[z]
               for ii = -1:IMAX[z]
                  u[z][ii, jj, kk, 1, c] = 1.0
                  u[z][ii, jj, kk, 2, c] = 0.0
                  u[z][ii, jj, kk, 3, c] = 0.0
                  u[z][ii, jj, kk, 4, c] = 0.0
                  u[z][ii, jj, kk, 5, c] = 1.0
               end
            end
         end
      end

#---------------------------------------------------------------------
# first store the "interpolated" values everywhere on the grid    
#---------------------------------------------------------------------
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

                   for ix = 1:2
                     Pface1[ix] = exact_solution(float(ix-1), eta, zeta)
                   end

                   for iy = 1:2
                     Pface2[iy] = exact_solution(xi, float(iy-1) , zeta)
                   end
 
                   for iz = 1:2
                     Pface3[iz] = exact_solution(xi, eta, float(iz-1))
                    # @info "$(xi), $(eta), $(float(iz-1)) --- $(Pface3[iz])"
                   end

                   for m = 1:5
                      Pxi   = xi   * Pface1[2][m] +(1.0e0-xi)   * Pface1[1][m]
                      Peta  = eta  * Pface2[2][m] +(1.0e0-eta)  * Pface2[1][m]
                      Pzeta = zeta * Pface3[2][m] +(1.0e0-zeta) * Pface3[1][m]

                      u[z][ii, jj, kk, m, c] = Pxi + Peta + Pzeta - Pxi*Peta - Pxi*Pzeta - Peta*Pzeta + Pxi*Peta*Pzeta
                   end
                   ii = ii + 1
                end
                jj = jj + 1
             end
             kk = kk+1
          end
       end

#---------------------------------------------------------------------
# now store the exact values on the boundaries        
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# west face                                                  
#---------------------------------------------------------------------
       c = slice[z][1, 1]
       ii = 0
       xi = 0.0e0
       kk = 0
       for k = cell_low[z][3, c]:cell_high[z][3, c]
          zeta = float(k) * dnzm1
          jj = 0
          for j = cell_low[z][2, c]:cell_high[z][2, c]
             eta = float(j) * dnym1
             temp = exact_solution(xi, eta, zeta)
             for m = 1:5
                u[z][ii, jj, kk, m, c] = temp[m]
             end
             jj = jj + 1
          end
          kk = kk + 1
       end

#---------------------------------------------------------------------
# east face                                                      
#---------------------------------------------------------------------
       c  = slice[z][1, ncells]
       ii = cell_size[z][1, c]-1
       xi = 1.0e0
       kk = 0
       for k = cell_low[z][3, c]:cell_high[z][3, c]
          zeta = float(k) * dnzm1
          jj = 0
          for j = cell_low[z][2, c]:cell_high[z][2, c]
             eta = float(j) * dnym1
             temp = exact_solution(xi, eta, zeta)
             for m = 1:5
                u[z][ii, jj, kk, m, c] = temp[m]
             end
             jj = jj + 1
          end
          kk = kk + 1
       end

#---------------------------------------------------------------------
# south face                                                 
#---------------------------------------------------------------------
       c = slice[z][2, 1]
       jj = 0
       eta = 0.0e0
       kk = 0
       for k = cell_low[z][3, c]:cell_high[z][3, c]
          zeta = float(k) * dnzm1
          ii = 0
          for i = cell_low[z][1, c]:cell_high[z][1, c]
             xi = float(i) * dnxm1
             temp = exact_solution(xi, eta, zeta)
             for m = 1:5
                u[z][ii, jj, kk, m, c] = temp[m]
             end
             ii = ii + 1
          end
          kk = kk + 1
       end


#---------------------------------------------------------------------
# north face                                    
#---------------------------------------------------------------------
       c = slice[z][2, ncells]
       jj = cell_size[z][2, c]-1
       eta = 1.0e0
       kk = 0
       for k = cell_low[z][3, c]:cell_high[z][3, c]
          zeta = float(k) * dnzm1
          ii = 0
          for i = cell_low[z][1, c]:cell_high[z][1, c]
             xi = float(i) * dnxm1
             temp = exact_solution(xi, eta, zeta)
             for m = 1:5
                u[z][ii, jj, kk, m, c] = temp[m]
             end
             ii = ii + 1
          end
          kk = kk + 1
       end

#---------------------------------------------------------------------
# bottom face                                       
#---------------------------------------------------------------------
       c = slice[z][3, 1]
       kk = 0
       zeta = 0.0e0
       jj = 0
       for j = cell_low[z][2, c]:cell_high[z][2, c]
          eta = float(j) * dnym1
          ii = 0
          for i = cell_low[z][1, c]:cell_high[z][1, c]
             xi = float(i) *dnxm1
             temp = exact_solution(xi, eta, zeta)
             for m = 1:5
                u[z][ii, jj, kk, m, c] = temp[m]
             end
             ii = ii + 1
          end
          jj = jj + 1
       end

#---------------------------------------------------------------------
# top face     
#---------------------------------------------------------------------
       c = slice[z][3, ncells]
       kk = cell_size[z][3, c]-1
       zeta = 1.0e0
       jj = 0
       for j = cell_low[z][2, c]:cell_high[z][2, c]
          eta = float(j) * dnym1
          ii = 0
          for i = cell_low[z][1, c]:cell_high[z][1, c]
             xi = float(i) * dnxm1
             temp = exact_solution(xi, eta, zeta)
             for m = 1:5
                u[z][ii, jj, kk, m, c] = temp[m]
             end
             ii = ii + 1
          end
          jj = jj + 1
       end

       return nothing
end


function lhsinit(z)

#---------------------------------------------------------------------
# loop over all cells                                       
#---------------------------------------------------------------------
       for c = 1:ncells

#---------------------------------------------------------------------
#         first, initialize the start and end arrays
#---------------------------------------------------------------------
          for d = 1:3
             if cell_coord[z][d, c] == 1
                cell_start[z][d, c] = 1
             else
                cell_start[z][d, c] = 0
             end
             if cell_coord[z][d, c] == ncells
                cell_end[z][d, c] = 1
             else
                cell_end[z][d, c] = 0
             end
          end

#---------------------------------------------------------------------
#     zap the whole left hand side for starters
#---------------------------------------------------------------------
          for n = 1:15
             for k = 0:cell_size[z][3, c]-1
                for j = 0:cell_size[z][2, c]-1
                   for i = 0:cell_size[z][1, c]-1
                      lhs[z][i, j, k, n, c] = 0.0e0
                   end
                end
             end
          end

#---------------------------------------------------------------------
# next, set all diagonal values to 1. This is overkill, but convenient
#---------------------------------------------------------------------
          for n = 1:3
             for k = 0:cell_size[z][3, c]-1
                for j = 0:cell_size[z][2, c]-1
                   for i = 0:cell_size[z][1, c]-1
                      lhs[z][i, j, k, 5*n-2, c] = 1.0e0
                   end
                end
             end
          end

       end

      return nothing
end



