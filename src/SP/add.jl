
#---------------------------------------------------------------------
# addition of update to the vector u
#---------------------------------------------------------------------

function add()

       for c = 1:ncells
          for m = 1:5
             for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
                for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
                   for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                      u[i, j, k, m, c] += rhs[i, j, k, m, c]
                   end
                end
             end
          end
       end 

 # VECTORIZATION ATTEMP (does not work)
 #=    for c = 1:ncells
         uu = view(u, cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1, 
                      cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1,
                      cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1,
                      1:5,
                      c)
         rhsrhs = view(rhs,cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1, 
                           cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1,
                           cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1,
                           1:5,
                           c)
         uu += rhsrhs
       end
=#

   

       return nothing
end
