#---------------------------------------------------------------------
#   set the boundary values of dependent variables
#---------------------------------------------------------------------

 function setbv()

#---------------------------------------------------------------------
#   set the dependent variable values along the top and bottom faces
#---------------------------------------------------------------------
      for j = 1:ny
         jglob = jpt + j
         for i = 1:nx
           iglob = ipt + i
            exact( iglob, jglob, 1, view(u, 1:5, i, j, 1 ) )
            exact( iglob, jglob, nz, view(u, 1:5, i, j, nz ) )
         end
      end

#---------------------------------------------------------------------
#   set the dependent variable values along north and south faces
#---------------------------------------------------------------------
      if west == -1
         for k = 1:nz
            for i = 1:nx
               iglob = ipt + i
               exact( iglob, 1, k, view(u, 1:5, i, 1, k ) )
            end
         end
      end

      if east == -1
          for k = 1:nz
             for i = 1:nx
                iglob = ipt + i
                exact( iglob, ny0, k, view(u, 1:5, i, ny, k ) )
             end
          end
      end

#---------------------------------------------------------------------
#   set the dependent variable values along east and west faces
#---------------------------------------------------------------------
      if north == -1
         for k = 1:nz
            for j = 1:ny
               jglob = jpt + j
               exact( 1, jglob, k, view(u, 1:5, 1, j, k))
            end
         end
      end

      if south == -1
         for k = 1:nz
            for j = 1:ny
               jglob = jpt + j
               exact( nx0, jglob, k, view(u, 1:5, nx, j, k) )
            end
         end
      end

      return nothing
end
