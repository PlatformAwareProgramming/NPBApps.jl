#---------------------------------------------------------------------
#     This function allocates space for a set of cells and fills the set     
#     such that communication between cells on different nodes is only
#     nearest neighbor                                                   
#---------------------------------------------------------------------

function make_set(node, no_nodes, grid_points, cell_coord, cell_low, cell_high, cell_size, slice)

#---------------------------------------------------------------------
#     compute square root; add small number to allow for roundoff
#     (note: this is computed in setup_mpi.f also, but prefer to do
#     it twice because of some include file problems).
#---------------------------------------------------------------------
      ncells = Int(sqrt(no_nodes))

#---------------------------------------------------------------------
#     this makes coding easier
#---------------------------------------------------------------------
      p = ncells

#---------------------------------------------------------------------
#     determine the location of the cell at the bottom of the 3D 
#     array of cells
#---------------------------------------------------------------------
      cell_coord[1, 1] = mod(node,p)
      cell_coord[2, 1] = div(node,p)
      cell_coord[3, 1] = 0

#---------------------------------------------------------------------
#     set the cell_coords for cells in the rest of the z-layers; 
#     this comes down to a simple linear numbering in the z-direct-
#     ion, and to the doubly-cyclic numbering in the other dirs     
#---------------------------------------------------------------------
      for c = 2:p
         cell_coord[1, c] = mod(cell_coord[1, c-1]+1, p)
         cell_coord[2, c] = mod(cell_coord[2, c-1]-1+p, p)
         cell_coord[3, c] = c-1
      end

#---------------------------------------------------------------------
#     offset all the coordinates by 1 to adjust for Fortran arrays
#---------------------------------------------------------------------
      for dir = 1:3
         for c = 1:p
            cell_coord[dir, c] = cell_coord[dir, c] + 1
         end
      end

#---------------------------------------------------------------------
#     slice[dir,n] contains the sequence number of the cell that is in
#     coordinate plane n in the dir direction
#---------------------------------------------------------------------
      for dir = 1:3
         for c = 1:p
            slice[dir, cell_coord[dir, c]] = c
         end
      end


#---------------------------------------------------------------------
#     fill the predecessor and successor entries, using the indices 
#     of the bottom cells (they are the same at each level of k 
#     anyway) acting as if full periodicity pertains; note that p is
#     added to those arguments to the mod functions that might
#     otherwise return wrong values when using the modulo function
#---------------------------------------------------------------------
      i = cell_coord[1, 1]-1
      j = cell_coord[2, 1]-1

      predecessor[1] = mod(i-1+p, p) + p*j
      predecessor[2] = i + p*mod(j-1+p, p)
      predecessor[3] = mod(i+1, p) + p*mod(j-1+p, p)
      successor[1]   = mod(i+1, p) + p*j
      successor[2]   = i + p*mod(j+1, p)
      successor[3]   = mod(i-1+p, p) + p*mod(j+1, p)

#---------------------------------------------------------------------
#     now compute the sizes of the cells                                    
#---------------------------------------------------------------------
      for dir = 1:3
#---------------------------------------------------------------------
#     set cell_coord range for each direction                            
#---------------------------------------------------------------------
         SIZE   = div(grid_points[dir],p)
         excess = mod(grid_points[dir], p)
         for c = 1:ncells
            if cell_coord[dir, c] <= excess
               cell_size[dir, c] = SIZE+1
               cell_low[dir, c] = (cell_coord[dir, c]-1)*(SIZE+1)
               cell_high[dir, c] = cell_low[dir, c]+SIZE
            else
               cell_size[dir, c] = SIZE
               cell_low[dir, c]  = excess*(SIZE+1)+(cell_coord[dir, c]-excess-1)*SIZE
               cell_high[dir, c] = cell_low[dir, c]+SIZE-1
            end
            if cell_size[dir, c] <= 2
               @printf(stdout, " Error: Cell size too small. Min size is 3\n", )
               ierrcode = 1
               MPI.Abort(MPI.COMM_WORLD, ierrcode)
               exit(1)
            end
         end
      end

      return ncells
end

#---------------------------------------------------------------------
#---------------------------------------------------------------------


