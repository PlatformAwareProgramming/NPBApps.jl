


#---------------------------------------------------------------------
#---------------------------------------------------------------------

 function neighbors(row, col, xdim, ydim)

#---------------------------------------------------------------------
#     figure out the neighbors and their wrap numbers for each processor
#---------------------------------------------------------------------

      south = -1
      east  = -1
      north = -1
      west  = -1

      if row > 1
          north = id -1
      else
          north = -1
      end

      if row < xdim
          south = id + 1
      else
          south = -1
      end

      if col > 1
          west = id- xdim
      else
          west = -1
      end

      if col < ydim
          east = id + xdim
      else
          east = -1
      end

      return south, east, north, west
end
