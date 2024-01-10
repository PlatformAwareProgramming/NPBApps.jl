#---------------------------------------------------------------------
#     addition of update to the vector u
#---------------------------------------------------------------------

function add(
            cell_size,
            cell_start,
            cell_end,         
            u,
            rhs,
            ::Val{ncells},   
            ) where ncells

      for c = 1:ncells
         for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
            for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
               for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                  for m = 1:5
                     u[m, i, j, k, c] += rhs[m, i, j, k, c]
                  end
               end
            end
         end
      end

       # VECTORIZATION ATTEMPT (does not work)
     #=for c = 1:ncells
         uu = view(u, 1:5,
                      cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1, 
                      cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1,
                      cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1,                      
                      c)
         rhsrhs = view(rhs, 1:5, 
                           cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1, 
                           cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1,
                           cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1,
                           c)
         uu .+= rhsrhs
       end=#


      return nothing
end
