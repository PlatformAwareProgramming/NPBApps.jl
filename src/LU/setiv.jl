using FortranFiles
using OffsetArrays
using Parameters
using Printf


#---------------------------------------------------------------------
#---------------------------------------------------------------------
      function setiv()

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#
#   set the initial values of independent variables based on tri-linear
#   interpolation of boundary values in the computational space.
#
#---------------------------------------------------------------------

#      use lu_data
#      implicit none

#---------------------------------------------------------------------
#  local variables
#---------------------------------------------------------------------
#      integer i, j, k, m
#      integer iglob, jglob
#      DOUBLEPRECISION  xi, eta, zeta
#      DOUBLEPRECISION  pxi, peta, pzeta
#      DOUBLEPRECISION  ue_1jk[5],ue_nx0jk[5],ue_i1k[5],  
#              ue_iny0k[5],ue_ij1[5],ue_ijnz[5]


      for k = 2:nz - 1
         zeta = ( float(k-1) ) / (nz-1)
         for j = 1:ny
          jglob = jpt + j
          if jglob != 1 && jglob != ny0
            eta = ( float(jglob-1) ) / (ny0-1)
            for i = 1:nx
              iglob = ipt + i
              if iglob != 1 && iglob != nx0
               xi = ( float(iglob-1) ) / (nx0-1)
               exact(1, jglob, k, ue_1jk)
               exact(nx0, jglob, k, ue_nx0jk)
               exact(iglob, 1, k, ue_i1k)
               exact(iglob, ny0, k, ue_iny0k)
               exact(iglob, jglob, 1, ue_ij1)
               exact(iglob, jglob, nz, ue_ijnz)
               for m = 1:5
                  pxi =   ( 1.0e+00 - xi ) * ue_1jk[m]+
                                     xi   * ue_nx0jk[m]
                  peta =  ( 1.0e+00 - eta ) * ue_i1k[m]+
                                     eta   * ue_iny0k[m]
                  pzeta = ( 1.0e+00 - zeta ) * ue_ij1[m]+
                                     zeta   * ue_ijnz[m]

                  u[ m, i, j, k ] = pxi + peta + pzeta-
                        pxi * peta - peta * pzeta - pzeta * pxi+
                        pxi * peta * pzeta

               end
              end
            end
          end
         end
      end

      return nothing
      end
