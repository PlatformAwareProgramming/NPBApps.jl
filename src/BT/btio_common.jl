

#---------------------------------------------------------------------
#---------------------------------------------------------------------

 function clear_timestep()

      for cio = 1:ncells
          for kio = 0:cell_size[z][3, cio]-1
              for jio = 0:cell_size[z][2, cio]-1
                  for ix = 0:cell_size[z][1, cio]-1
                            u[1, ix, jio, kio, cio] = 0
                            u[2, ix, jio, kio, cio] = 0
                            u[3, ix, jio, kio, cio] = 0
                            u[4, ix, jio, kio, cio] = 0
                            u[5, ix, jio, kio, cio] = 0
                  end
              end
          end
      end

      return nothing
end

