ue_1jk = Array{FloatType}(undef, 5)
ue_nx0jk = Array{FloatType}(undef, 5)
ue_i1k = Array{FloatType}(undef, 5) 
ue_iny0k = Array{FloatType}(undef, 5)
ue_ij1 = Array{FloatType}(undef, 5)
ue_ijnz = Array{FloatType}(undef, 5)

#---------------------------------------------------------------------
#
#   set the initial values of independent variables based on tri-linear
#   interpolation of boundary values in the computational space.
#
#---------------------------------------------------------------------

function setiv(u, nx0, ny0, nz0, nx, ny, nz, ipt, jpt)

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

function write_u2(z, u, nx0, ny0, nz0, nx, ny, nz, ipt, jpt)

  for k = 2:nz - 1
    for j = 1:ny
#     jglob = jpt + j
#     if jglob != 1 && jglob != ny0
       for i = 1:nx
#         iglob = ipt + i
#         if iglob != 1 && iglob != nx0
           for m = 1:5
              @info "$z, $m, $i, $j, $k, $(u[m, i, j, k])"
           end
          #end
       end
#      end
    end
 end

end

function write_u(z, u, nx0, ny0, nz0, nx, ny, nz, ipt, jpt)

  for k = 1:nz
    for j = 1:ny
     jglob = jpt + j
#     if jglob != 1 && jglob != ny0
       for i = 1:nx
         iglob = ipt + i
#         if iglob != 1 && iglob != nx0
           for m = 1:5
              @info "$node, $z, $m, $iglob, $jglob, $k, $(u[m, i, j, k])"
           end
          #end
       end
#      end
    end
 end

end
