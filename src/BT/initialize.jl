using FortranFiles
using OffsetArrays
using Parameters
using Printf

#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function initialize()

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     This subroutine initializes the field variable u using 
#     tri-linear transfinite interpolation of the boundary values     
#---------------------------------------------------------------------

#      use bt_data
#      implicit none

#      integer c, i, j, k, m, ii, jj, kk, ix, iy, iz
#      DOUBLEPRECISION  xi, eta, zeta, Pface[5,3,2], Pxi, Peta,  
#           Pzeta, temp[5]

#---------------------------------------------------------------------
#  Later (in compute_rhs) we compute 1/u for every element. A few of 
#  the corner elements are not used, but it convenient (and faster) 
#  to compute the whole thing with a simple loop. Make sure those 
#  values are nonzero by initializing the whole thing here. 
#---------------------------------------------------------------------
      for c = 1:ncells
         for kk = -1:KMAX
            for jj = -1:JMAX
               for ii = -1:IMAX
                  for m = 1:5
                     u[m, ii, jj, kk, c] = 1.0
                  end
               end
            end
         end
      end
#---------------------------------------------------------------------



#---------------------------------------------------------------------
#     first store the "interpolated" values everywhere on the grid    
#---------------------------------------------------------------------
      for c = 1:ncells
         kk = 0
         for k = cell_low[3, c]:cell_high[3, c]
            zeta = float(k) * dnzm1
            jj = 0
            for j = cell_low[2, c]:cell_high[2, c]
               eta = float(j) * dnym1
               ii = 0
               for i = cell_low[1, c]:cell_high[1, c]
                  xi = float(i) * dnxm1

                  for ix = 1:2
                     exact_solution(float(ix-1), eta, zeta,
                          Pface[1,1,ix])
                  end

                  for iy = 1:2
                     exact_solution(xi, float(iy-1) , zeta,
                          Pface[1,2,iy])
                  end

                  for iz = 1:2
                     exact_solution(xi, eta, float(iz-1),
                          Pface[1,3,iz])
                  end

                  for m = 1:5
                     Pxi   = xi   * Pface[m,1,2] +(
                          1.0e0-xi)   * Pface[m,1,1]
                     Peta  = eta  * Pface[m,2,2] +(
                          1.0e0-eta)  * Pface[m,2,1]
                     Pzeta = zeta * Pface[m,3,2] +(
                          1.0e0-zeta) * Pface[m,3,1]

                     u[m, ii, jj, kk, c] = Pxi + Peta + Pzeta -
                          Pxi*Peta - Pxi*Pzeta - Peta*Pzeta +
                          Pxi*Peta*Pzeta

                  end
                  ii = ii + 1
               end
               jj = jj + 1
            end
            kk = kk+1
         end
      end

#---------------------------------------------------------------------
#     now store the exact values on the boundaries        
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#     west face                                                  
#---------------------------------------------------------------------
      c = slice[1, 1]
      ii = 0
      xi = 0.0e0
      kk = 0
      for k = cell_low[3, c]:cell_high[3, c]
         zeta = float(k) * dnzm1
         jj = 0
         for j = cell_low[2, c]:cell_high[2, c]
            eta = float(j) * dnym1
            exact_solution(xi, eta, zeta, temp)
            for m = 1:5
               u[m, ii, jj, kk, c] = temp[m]
            end
            jj = jj + 1
         end
         kk = kk + 1
      end

#---------------------------------------------------------------------
#     east face                                                      
#---------------------------------------------------------------------
      c  = slice[1, ncells]
      ii = cell_size[1, c]-1
      xi = 1.0e0
      kk = 0
      for k = cell_low[3, c]:cell_high[3, c]
         zeta = float(k) * dnzm1
         jj = 0
         for j = cell_low[2, c]:cell_high[2, c]
            eta = float(j) * dnym1
            exact_solution(xi, eta, zeta, temp)
            for m = 1:5
               u[m, ii, jj, kk, c] = temp[m]
            end
            jj = jj + 1
         end
         kk = kk + 1
      end

#---------------------------------------------------------------------
#     south face                                                 
#---------------------------------------------------------------------
      c = slice[2, 1]
      jj = 0
      eta = 0.0e0
      kk = 0
      for k = cell_low[3, c]:cell_high[3, c]
         zeta = float(k) * dnzm1
         ii = 0
         for i = cell_low[1, c]:cell_high[1, c]
            xi = float(i) * dnxm1
            exact_solution(xi, eta, zeta, temp)
            for m = 1:5
               u[m, ii, jj, kk, c] = temp[m]
            end
            ii = ii + 1
         end
         kk = kk + 1
      end


#---------------------------------------------------------------------
#     north face                                    
#---------------------------------------------------------------------
      c = slice[2, ncells]
      jj = cell_size[2, c]-1
      eta = 1.0e0
      kk = 0
      for k = cell_low[3, c]:cell_high[3, c]
         zeta = float(k) * dnzm1
         ii = 0
         for i = cell_low[1, c]:cell_high[1, c]
            xi = float(i) * dnxm1
            exact_solution(xi, eta, zeta, temp)
            for m = 1:5
               u[m, ii, jj, kk, c] = temp[m]
            end
            ii = ii + 1
         end
         kk = kk + 1
      end

#---------------------------------------------------------------------
#     bottom face                                       
#---------------------------------------------------------------------
      c = slice[3, 1]
      kk = 0
      zeta = 0.0e0
      jj = 0
      for j = cell_low[2, c]:cell_high[2, c]
         eta = float(j) * dnym1
         ii = 0
         for i = cell_low[1, c]:cell_high[1, c]
            xi = float(i) *dnxm1
            exact_solution(xi, eta, zeta, temp)
            for m = 1:5
               u[m, ii, jj, kk, c] = temp[m]
            end
            ii = ii + 1
         end
         jj = jj + 1
      end

#---------------------------------------------------------------------
#     top face     
#---------------------------------------------------------------------
      c = slice[3, ncells]
      kk = cell_size[3, c]-1
      zeta = 1.0e0
      jj = 0
      for j = cell_low[2, c]:cell_high[2, c]
         eta = float(j) * dnym1
         ii = 0
         for i = cell_low[1, c]:cell_high[1, c]
            xi = float(i) * dnxm1
            exact_solution(xi, eta, zeta, temp)
            for m = 1:5
               u[m, ii, jj, kk, c] = temp[m]
            end
            ii = ii + 1
         end
         jj = jj + 1
      end

      return nothing
      end


#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function lhsinit()

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#      use bt_data
#      implicit none

#      integer i, j, k, d, c, m, n

#---------------------------------------------------------------------
#     loop over all cells                                       
#---------------------------------------------------------------------
      for c = 1:ncells

#---------------------------------------------------------------------
#     first, initialize the start and end arrays
#---------------------------------------------------------------------
         for d = 1:3
            if cell_coord[d, c] == 1
               cell_start[d, c] = 1
            else
               cell_start[d, c] = 0
            end
            if cell_coord[d, c] == ncells
               cell_end[d, c] = 1
            else
               cell_end[d, c] = 0
            end
         end

#---------------------------------------------------------------------
#     zero the whole left hand side for starters
#---------------------------------------------------------------------
         for k = 0:cell_size[3, c]-1
            for j = 0:cell_size[2, c]-1
               for i = 0:cell_size[1, c]-1
                  for m = 1:5
                     for n = 1:5
                         lhsc[m, n, i, j, k, c] = 0.0e0
                     end
                  end
               end
            end
         end

      end

      return nothing
      end


#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function lhsabinit(lhsa, lhsb, SIZE)
#      implicit none

#      integer SIZE
#      DOUBLEPRECISION lhsa[5, 5, -1:SIZE], lhsb[5, 5, -1:SIZE]

#      integer i, m, n

#---------------------------------------------------------------------
#     next, set all diagonal values to 1. This is overkill, but convenient
#---------------------------------------------------------------------
      for i = 0:SIZE
         for m = 1:5
            for n = 1:5
               lhsa[m,n,i] = 0.0e0
               lhsb[m,n,i] = 0.0e0
            end
            lhsb[m,m,i] = 1.0e0
         end
      end

      return nothing
      end



