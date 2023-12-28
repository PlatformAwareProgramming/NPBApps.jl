#---------------------------------------------------------------------
# This function computes the left hand side for the three y-factors   
#---------------------------------------------------------------------

function lhsy(c,
            cell_size,
            cell_start,
            cell_end,
            lhs,
            rho_i,
            vs,
            rhoq,
            cv,
            speed,
            con43, c3c4, c1c5, c2dtty1,
            dy3, dy5, dy1, dymax, dtty1, dtty2, 
            comz5, comz4, comz1, comz6)

#---------------------------------------------------------------------
#      treat only cell c
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#      first fill the lhs for the u-eigenvalue         
#---------------------------------------------------------------------
       for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
          for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1

             for j = cell_start[2, c]-1:cell_size[2, c]-cell_end[2, c]
                ru1 = c3c4*rho_i[i, j, k, c]
                cv[j] = vs[i, j, k, c]
                rhoq[j] = max( dy3 + con43 * ru1,
                                 dy5 + c1c5*ru1,
                                 dymax + ru1,
                                 dy1)
             end

             for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
                lhs[i, j, k, 1, c] =  0.0e0
                lhs[i, j, k, 2, c] = -dtty2 * cv[j-1] - dtty1 * rhoq[j-1]
                lhs[i, j, k, 3, c] =  1.0 + c2dtty1 * rhoq[j]
                lhs[i, j, k, 4, c] =  dtty2 * cv[j+1] - dtty1 * rhoq[j+1]
                lhs[i, j, k, 5, c] =  0.0e0
             end
          end
       end

#---------------------------------------------------------------------
#      add fourth order dissipation                             
#---------------------------------------------------------------------
       if cell_start[2, c] > 0
          j = 1
          for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
             for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1

                lhs[i, j, k, 3, c] += comz5
                lhs[i, j, k, 4, c] -= comz4
                lhs[i, j, k, 5, c] += comz1

                lhs[i, j+1, k, 2, c] -= comz4
                lhs[i, j+1, k, 3, c] += comz6
                lhs[i, j+1, k, 4, c] -= comz4
                lhs[i, j+1, k, 5, c] += comz1
             end
          end
       end

       for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
          for j = 3*cell_start[2, c]:cell_size[2, c]-3*cell_end[2, c]-1
             for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1

                lhs[i, j, k, 1, c] += comz1
                lhs[i, j, k, 2, c] -= comz4
                lhs[i, j, k, 3, c] += comz6
                lhs[i, j, k, 4, c] -= comz4
                lhs[i, j, k, 5, c] += comz1
             end
          end
       end

       if cell_end[2, c] > 0
          j = cell_size[2, c]-3
          for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
             for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                lhs[i, j, k, 1, c] += comz1
                lhs[i, j, k, 2, c] -= comz4
                lhs[i, j, k, 3, c] += comz6
                lhs[i, j, k, 4, c] -= comz4

                lhs[i, j+1, k, 1, c] += comz1
                lhs[i, j+1, k, 2, c] -= comz4
                lhs[i, j+1, k, 3, c] += comz5
             end
          end
       end

#---------------------------------------------------------------------
#      subsequently, do the other two factors                    
#---------------------------------------------------------------------
       for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
          for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
             for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                lhs[i, j, k, 1+5, c]  = lhs[i, j, k, 1, c]
                lhs[i, j, k, 2+5, c]  = lhs[i, j, k, 2, c] - dtty2 * speed[i, j-1, k, c]
                lhs[i, j, k, 3+5, c]  = lhs[i, j, k, 3, c]
                lhs[i, j, k, 4+5, c]  = lhs[i, j, k, 4, c] + dtty2 * speed[i, j+1, k, c]
                lhs[i, j, k, 5+5, c] = lhs[i, j, k, 5, c]
                lhs[i, j, k, 1+10, c] = lhs[i, j, k, 1, c]
                lhs[i, j, k, 2+10, c] = lhs[i, j, k, 2, c] + dtty2 * speed[i, j-1, k, c]
                lhs[i, j, k, 3+10, c] = lhs[i, j, k, 3, c]
                lhs[i, j, k, 4+10, c] = lhs[i, j, k, 4, c] - dtty2 * speed[i, j+1, k, c]
                lhs[i, j, k, 5+10, c] = lhs[i, j, k, 5, c]
             end
          end
       end

       return nothing
end



