using FortranFiles
using OffsetArrays
using Parameters
using Printf

#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function add()

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     addition of update to the vector u
#---------------------------------------------------------------------

#      use bt_data
#      implicit none

#      integer  c, i, j, k, m

      for c = 1:ncells
         for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
            for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
               for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                  for m = 1:5
                     u[m, i, j, k, c] = u[m, i, j, k, c] + rhs[m, i, j, k, c]
                  end
               end
            end
         end
      end

      return nothing
      end
