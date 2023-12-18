#---------------------------------------------------------------------
# This function computes the left hand side for the three x-factors  
#---------------------------------------------------------------------

function lhsx(c)

#---------------------------------------------------------------------
#      treat only cell c             
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#      first fill the lhs for the u-eigenvalue                   
#---------------------------------------------------------------------
       for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
          for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
             for i = cell_start[1, c]-1:cell_size[1, c]-cell_end[1, c]
                ru1 = c3c4*rho_i[i, j, k, c]
                cv[i] = us[i, j, k, c]
                rhon[i] = max(dx2+con43*ru1, dx5+c1c5*ru1, dxmax+ru1, dx1)
             end

             for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                lhs[i, j, k, 1, c] =   0.0e0
                lhs[i, j, k, 2, c] = - dttx2 * cv[i-1] - dttx1 * rhon[i-1]
                lhs[i, j, k, 3, c] =   1.0e0 + c2dttx1 * rhon[i]
                lhs[i, j, k, 4, c] =   dttx2 * cv[i+1] - dttx1 * rhon[i+1]
                lhs[i, j, k, 5, c] =   0.0e0
             end
          end
       end

#---------------------------------------------------------------------
#      add fourth order dissipation                             
#---------------------------------------------------------------------
       if cell_start[1, c] > 0
          i = 1
          for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
             for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
                lhs[i, j, k, 3, c] = lhs[i, j, k, 3, c] + comz5
                lhs[i, j, k, 4, c] = lhs[i, j, k, 4, c] - comz4
                lhs[i, j, k, 5, c] = lhs[i, j, k, 5, c] + comz1

                lhs[i+1, j, k, 2, c] = lhs[i+1, j, k, 2, c] - comz4
                lhs[i+1, j, k, 3, c] = lhs[i+1, j, k, 3, c] + comz6
                lhs[i+1, j, k, 4, c] = lhs[i+1, j, k, 4, c] - comz4
                lhs[i+1, j, k, 5, c] = lhs[i+1, j, k, 5, c] + comz1
             end
          end
       end

       for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
          for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
             for i = 3*cell_start[1, c]:cell_size[1, c]-3*cell_end[1, c]-1
                lhs[i, j, k, 1, c] = lhs[i, j, k, 1, c] + comz1
                lhs[i, j, k, 2, c] = lhs[i, j, k, 2, c] - comz4
                lhs[i, j, k, 3, c] = lhs[i, j, k, 3, c] + comz6
                lhs[i, j, k, 4, c] = lhs[i, j, k, 4, c] - comz4
                lhs[i, j, k, 5, c] = lhs[i, j, k, 5, c] + comz1
             end
          end
       end

       if cell_end[1, c] > 0
          i = cell_size[1, c]-3
          for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
             for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
                lhs[i, j, k, 1, c] = lhs[i, j, k, 1, c] + comz1
                lhs[i, j, k, 2, c] = lhs[i, j, k, 2, c] - comz4
                lhs[i, j, k, 3, c] = lhs[i, j, k, 3, c] + comz6
                lhs[i, j, k, 4, c] = lhs[i, j, k, 4, c] - comz4

                lhs[i+1, j, k, 1, c] = lhs[i+1, j, k, 1, c] + comz1
                lhs[i+1, j, k, 2, c] = lhs[i+1, j, k, 2, c] - comz4
                lhs[i+1, j, k, 3, c] = lhs[i+1, j, k, 3, c] + comz5
             end
          end
       end

#---------------------------------------------------------------------
#      subsequently, fill the other factors (u+c), (u-c) by a4ing to 
#      the first  
#---------------------------------------------------------------------
       for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
          for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
             for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                lhs[i, j, k, 1+5, c]  = lhs[i, j, k, 1, c]
                lhs[i, j, k, 2+5, c]  = lhs[i, j, k, 2, c] -
                                  dttx2 * speed[i-1, j, k, c]
                lhs[i, j, k, 3+5, c]  = lhs[i, j, k, 3, c]
                lhs[i, j, k, 4+5, c]  = lhs[i, j, k, 4, c] +
                                  dttx2 * speed[i+1, j, k, c]
                lhs[i, j, k, 5+5, c] = lhs[i, j, k, 5, c]
                lhs[i, j, k, 1+10, c] = lhs[i, j, k, 1, c]
                lhs[i, j, k, 2+10, c] = lhs[i, j, k, 2, c] +
                                  dttx2 * speed[i-1, j, k, c]
                lhs[i, j, k, 3+10, c] = lhs[i, j, k, 3, c]
                lhs[i, j, k, 4+10, c] = lhs[i, j, k, 4, c] -
                                  dttx2 * speed[i+1, j, k, c]
                lhs[i, j, k, 5+10, c] = lhs[i, j, k, 5, c]
             end
          end
       end

       return nothing
end



