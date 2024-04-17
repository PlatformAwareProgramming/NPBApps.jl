


#---------------------------------------------------------------------
#---------------------------------------------------------------------

 function neighbors(row, col)

#---------------------------------------------------------------------
#     figure out the neighbors and their wrap numbers for each processor
#---------------------------------------------------------------------

      south = -1
      east  = -1
      north = -1
      west  = -1

      if row > 1
          north2 = north = id - 1
      else
          north = -1
          north2 = col*xdim - 1
      end

      if row < xdim
          south2 = south = id + 1
      else
          south = -1
          south2 = (col-1)*xdim
      end

      if col > 1
          west2 = west = id - xdim
      else
          west = -1
          west2 = (ydim-1)*ydim + row - 1    # 3*4 + 2 - 1 = 13
      end

      if col < ydim
          east2 = east = id + xdim
      else
          east = -1
          east2 = row - 1
      end

      return south, east, north, west, south2, east2, north2, west2
end
