#---------------------------------------------------------------------
#---------------------------------------------------------------------

       function txinvr()

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# block-diagonal matrix-vector multiplication                  
#---------------------------------------------------------------------

#       use sp_data
#       implicit none

#       integer c, i, j, k
#       DOUBLEPRECISION t1, t2, t3, ac, ru1, uu, vv, ww, r1, r2, r3,  
#                        r4, r5, ac2inv

#---------------------------------------------------------------------
#      loop over all cells owned by this node          
#---------------------------------------------------------------------
       for c = 1:ncells
          for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
             for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
                for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1

                   ru1 = rho_i[i, j, k, c]
                   uu = us[i, j, k, c]
                   vv = vs[i, j, k, c]
                   ww = ws[i, j, k, c]
                   ac = speed[i, j, k, c]
                   ac2inv = ainv[i, j, k, c]*ainv[i, j, k, c]

                   r1 = rhs[i, j, k, 1, c]
                   r2 = rhs[i, j, k, 2, c]
                   r3 = rhs[i, j, k, 3, c]
                   r4 = rhs[i, j, k, 4, c]
                   r5 = rhs[i, j, k, 5, c]

                   t1 = c2 * ac2inv * ( qs[i, j, k, c]*r1 - uu*r2  -
                        vv*r3 - ww*r4 + r5 )
                   t2 = bt * ru1 * ( uu * r1 - r2 )
                   t3 = ( bt * ru1 * ac ) * t1

                   rhs[i, j, k, 1, c] = r1 - t1
                   rhs[i, j, k, 2, c] = - ru1 * ( ww*r1 - r4 )
                   rhs[i, j, k, 3, c] =   ru1 * ( vv*r1 - r3 )
                   rhs[i, j, k, 4, c] = - t2 + t3
                   rhs[i, j, k, 5, c] =   t2 + t3
                end
             end
          end
       end

       return nothing
       end


