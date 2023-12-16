#---------------------------------------------------------------------
#---------------------------------------------------------------------

       function pinvr(c)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#   block-diagonal matrix-vector multiplication                       
#---------------------------------------------------------------------

#       use sp_data
#       implicit none

#       integer i, j, k, c
#       DOUBLEPRECISION r1, r2, r3, r4, r5, t1, t2

#---------------------------------------------------------------------
#      treat only one cell                                   
#---------------------------------------------------------------------
       for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
          for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
             for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1

                r1 = rhs[i, j, k, 1, c]
                r2 = rhs[i, j, k, 2, c]
                r3 = rhs[i, j, k, 3, c]
                r4 = rhs[i, j, k, 4, c]
                r5 = rhs[i, j, k, 5, c]

                t1 = bt * r1
                t2 = 0.5e0 * ( r4 + r5 )

                rhs[i, j, k, 1, c] =  bt * ( r4 - r5 )
                rhs[i, j, k, 2, c] = -r3
                rhs[i, j, k, 3, c] =  r2
                rhs[i, j, k, 4, c] = -t1 + t2
                rhs[i, j, k, 5, c] =  t1 + t2
             end
          end
       end

       return nothing
       end



