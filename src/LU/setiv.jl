ue_1jk = Array{Float64}(undef, 5)
ue_nx0jk = Array{Float64}(undef, 5)
ue_i1k = Array{Float64}(undef, 5) 
ue_iny0k = Array{Float64}(undef, 5)
ue_ij1 = Array{Float64}(undef, 5)
ue_ijnz = Array{Float64}(undef, 5)

#---------------------------------------------------------------------
#
#   set the initial values of independent variables based on tri-linear
#   interpolation of boundary values in the computational space.
#
#---------------------------------------------------------------------

function setiv(nx0, ny0, nz0)

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
               exact(1, jglob, k, ue_1jk, nx0, ny0, nz0)
               exact(nx0, jglob, k, ue_nx0jk, nx0, ny0, nz0)
               exact(iglob, 1, k, ue_i1k, nx0, ny0, nz0)
               exact(iglob, ny0, k, ue_iny0k, nx0, ny0, nz0)
               exact(iglob, jglob, 1, ue_ij1, nx0, ny0, nz0)
               exact(iglob, jglob, nz, ue_ijnz, nx0, ny0, nz0)
               for m = 1:5
                  pxi =   (1.0e+00 - xi) * ue_1jk[m] + xi * ue_nx0jk[m]
                  peta =  (1.0e+00 - eta) * ue_i1k[m] + eta * ue_iny0k[m]
                  pzeta = (1.0e+00 - zeta) * ue_ij1[m] + zeta * ue_ijnz[m]

                  u[ m, i, j, k ] = pxi + peta + pzeta - pxi * peta - peta * pzeta - pzeta * pxi + pxi * peta * pzeta
               end
              end
            end
          end
         end
      end

      return nothing
end
