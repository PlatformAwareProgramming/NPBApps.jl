
#---------------------------------------------------------------------
#---------------------------------------------------------------------

function compute_rhs(ncells,
                     cell_size,
                     cell_start,
                     cell_end,
                     u,
                     rhs,
                     rho_i,
                     us,
                     vs,
                     ws,
                     square,
                     qs,
                     ainv,
                     speed,
                     forcing,
                     dt,
                     tx2,
                     ty2,
                     tz2,
                     c1,
                     c2,
                     c1c2,
                     dx1tx1,
                     dx2tx1,
                     dx3tx1, 
                     dx4tx1,
                     dx5tx1,
                     dy1ty1,
                     dy2ty1,
                     dy3ty1,
                     dy4ty1,
                     dy5ty1,
                     dz1tz1,
                     dz2tz1,
                     dz3tz1,
                     dz4tz1,
                     dz5tz1,
                     xxcon2,
                     xxcon3,
                     xxcon4,
                     xxcon5,
                     yycon2,
                     yycon3,
                     yycon4,
                     yycon5,
                     zzcon2,
                     zzcon3,
                     zzcon4,
                     zzcon5,
                     dssp,
                     con43,
                     )

       if (timeron) timer_start(t_rhs) end
#---------------------------------------------------------------------
# loop over all cells owned by this node                           
#---------------------------------------------------------------------
       for c = 1:ncells
         
#---------------------------------------------------------------------
#         compute the reciprocal of density, and the kinetic energy, 
#         and the speed of sound. 
#---------------------------------------------------------------------

          for k = -1:cell_size[3, c]
             for j = -1:cell_size[2, c]
                for i = -1:cell_size[1, c]
                   rho_inv = 1.0e0/u[i, j, k, 1, c]
                   rho_i[i, j, k, c] = rho_inv
                   us[i, j, k, c] = u[i, j, k, 2, c] * rho_inv
                   vs[i, j, k, c] = u[i, j, k, 3, c] * rho_inv
                   ws[i, j, k, c] = u[i, j, k, 4, c] * rho_inv
                   square[i, j, k, c]     = 0.5e0* (
                              u[i, j, k, 2, c]*u[i, j, k, 2, c] +
                              u[i, j, k, 3, c]*u[i, j, k, 3, c] +
                              u[i, j, k, 4, c]*u[i, j, k, 4, c] ) * rho_inv
                   qs[i, j, k, c] = square[i, j, k, c] * rho_inv
#---------------------------------------------------------------------
#                  (don't need speed and ainx until the lhs computation)
#---------------------------------------------------------------------
                   aux = c1c2*rho_inv* (u[i, j, k, 5, c] - square[i, j, k, c])
                   aux = sqrt(aux)
                   speed[i, j, k, c] = aux
                   ainv[i, j, k, c]  = 1.0e0/aux
                end
             end
          end

#---------------------------------------------------------------------
# copy the exact forcing term to the right hand side;  because 
# this forcing term is known, we can store it on the whole of every 
# cell,  including the boundary                   
#---------------------------------------------------------------------

          for m = 1:5
             for k = 0:cell_size[3, c]-1
                for j = 0:cell_size[2, c]-1
                   for i = 0:cell_size[1, c]-1
                      rhs[i, j, k, m, c] = forcing[i, j, k, m, c]
                   end
                end
             end
          end

#---------------------------------------------------------------------
#         compute xi-direction fluxes 
#---------------------------------------------------------------------
          for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
             for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
                for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                   uijk = us[i, j, k, c]
                   up1  = us[i+1, j, k, c]
                   um1  = us[i-1, j, k, c]

                   rhs[i, j, k, 1, c] += dx1tx1 *(
                          u[i+1, j, k, 1, c] - 2.0e0*u[i, j, k, 1, c] +
                           u[i-1, j, k, 1, c]) -
                          tx2 * (u[i+1, j, k, 2, c] - u[i-1, j, k, 2, c])

                   rhs[i, j, k, 2, c] += dx2tx1 *(
                          u[i+1, j, k, 2, c] - 2.0e0*u[i, j, k, 2, c] +
                           u[i-1, j, k, 2, c]) +
                          xxcon2*con43 * (up1 - 2.0e0*uijk + um1) -
                          tx2 * (u[i+1, j, k, 2, c]*up1 -
                                 u[i-1, j, k, 2, c]*um1 +(
                                 u[i+1, j, k, 5, c]- square[i+1, j, k, c]-
                                  u[i-1, j, k, 5, c]+ square[i-1, j, k, c])*
                                  c2)

                   rhs[i, j, k, 3, c] += dx3tx1 *(
                          u[i+1, j, k, 3, c] - 2.0e0*u[i, j, k, 3, c] +
                           u[i-1, j, k, 3, c]) +
                          xxcon2 * (vs[i+1, j, k, c] - 2.0e0*vs[i, j, k, c] +
                                    vs[i-1, j, k, c]) -
                          tx2 * (u[i+1, j, k, 3, c]*up1 -
                                 u[i-1, j, k, 3, c]*um1)

                   rhs[i, j, k, 4, c] += dx4tx1 *(
                          u[i+1, j, k, 4, c] - 2.0e0*u[i, j, k, 4, c] +
                           u[i-1, j, k, 4, c]) +
                          xxcon2 * (ws[i+1, j, k, c] - 2.0e0*ws[i, j, k, c] +
                                    ws[i-1, j, k, c]) -
                          tx2 * (u[i+1, j, k, 4, c]*up1 -
                                 u[i-1, j, k, 4, c]*um1)

                   rhs[i, j, k, 5, c] += dx5tx1 *(
                          u[i+1, j, k, 5, c] - 2.0e0*u[i, j, k, 5, c] +
                           u[i-1, j, k, 5, c]) +
                          xxcon3 * (qs[i+1, j, k, c] - 2.0e0*qs[i, j, k, c] +
                                    qs[i-1, j, k, c]) +
                          xxcon4 * (up1*up1 -       2.0e0*uijk*uijk +
                                    um1*um1) +
                          xxcon5 * (u[i+1, j, k, 5, c]*rho_i[i+1, j, k, c] -
                                    2.0e0*u[i, j, k, 5, c]*rho_i[i, j, k, c] +
                                    u[i-1, j, k, 5, c]*rho_i[i-1, j, k, c]) -
                          tx2 * ( (c1*u[i+1, j, k, 5, c] -
                                   c2*square[i+1, j, k, c])*up1 -(
                                  c1*u[i-1, j, k, 5, c] -
                                   c2*square[i-1, j, k, c])*um1 )
                end
             end
          end

#---------------------------------------------------------------------
#         add fourth order xi-direction dissipation               
#---------------------------------------------------------------------
          if cell_start[1, c] > 0
             i = 1
             for m = 1:5
                for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
                   for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
                      rhs[i, j, k, m, c] -= dssp *(5.0e0*u[i, j, k, m, c] - 4.0e0*u[i+1, j, k, m, c] + u[i+2, j, k, m, c])
                   end
                end
             end

             i = 2
             for m = 1:5
                for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
                   for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
                      rhs[i, j, k, m, c] -= dssp *(-4.0e0*u[i-1, j, k, m, c] + 6.0e0*u[i, j, k, m, c] - 4.0e0*u[i+1, j, k, m, c] + u[i+2, j, k, m, c])
                   end
                end
             end
          end

          for m = 1:5
             for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
                for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
                   for i = 3*cell_start[1, c]:cell_size[1, c]-3*cell_end[1, c]-1
                      rhs[i, j, k, m, c] -= dssp *(u[i-2, j, k, m, c] - 4.0e0*u[i-1, j, k, m, c] + 6.0*u[i, j, k, m, c] - 4.0e0*u[i+1, j, k, m, c] + u[i+2, j, k, m, c] )
                      #if i==0
                        #@info "i=-2 !!! rhs[$z][$i, $j, $k, $m, $c] =  $(rhs[i, j, k, m, c]) -- $(u[i-2, j, k, m, c])"
                      #end
                   end
                end
             end
          end


          if cell_end[1, c] > 0
             i = cell_size[1, c]-3
             for m = 1:5
                for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
                   for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
                      rhs[i, j, k, m, c] -= dssp *(u[i-2, j, k, m, c] - 4.0e0*u[i-1, j, k, m, c] + 6.0e0*u[i, j, k, m, c] - 4.0e0*u[i+1, j, k, m, c] )
                   end
                end
             end

             i = cell_size[1, c]-2
             for m = 1:5
                for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
                   for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
                      rhs[i, j, k, m, c] -= dssp *(u[i-2, j, k, m, c] - 4.0e0*u[i-1, j, k, m, c] + 5.0e0*u[i, j, k, m, c] )
                   end
                end
             end
          end

#---------------------------------------------------------------------
#         compute eta-direction fluxes 
#---------------------------------------------------------------------
          for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
             for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
                for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                   vijk = vs[i, j, k, c]
                   vp1  = vs[i, j+1, k, c]
                   vm1  = vs[i, j-1, k, c]
                   rhs[i, j, k, 1, c] += dy1ty1 *(
                         u[i, j+1, k, 1, c] - 2.0e0*u[i, j, k, 1, c] +
                          u[i, j-1, k, 1, c]) -
                         ty2 * (u[i, j+1, k, 3, c] - u[i, j-1, k, 3, c])
                   rhs[i, j, k, 2, c] += dy2ty1 *(
                         u[i, j+1, k, 2, c] - 2.0e0*u[i, j, k, 2, c] +
                          u[i, j-1, k, 2, c]) +
                         yycon2 * (us[i, j+1, k, c] - 2.0e0*us[i, j, k, c] +
                                   us[i, j-1, k, c]) -
                         ty2 * (u[i, j+1, k, 2, c]*vp1 -
                                u[i, j-1, k, 2, c]*vm1)
                   rhs[i, j, k, 3, c] += dy3ty1 *(
                         u[i, j+1, k, 3, c] - 2.0e0*u[i, j, k, 3, c] +
                          u[i, j-1, k, 3, c]) +
                         yycon2*con43 * (vp1 - 2.0e0*vijk + vm1) -
                         ty2 * (u[i, j+1, k, 3, c]*vp1 -
                                u[i, j-1, k, 3, c]*vm1 +(
                                u[i, j+1, k, 5, c] - square[i, j+1, k, c] -
                                 u[i, j-1, k, 5, c] + square[i, j-1, k, c])*
                                c2)
                   rhs[i, j, k, 4, c] += dy4ty1 *(
                         u[i, j+1, k, 4, c] - 2.0e0*u[i, j, k, 4, c] +
                          u[i, j-1, k, 4, c]) +
                         yycon2 * (ws[i, j+1, k, c] - 2.0e0*ws[i, j, k, c] +
                                   ws[i, j-1, k, c]) -
                         ty2 * (u[i, j+1, k, 4, c]*vp1 -
                                u[i, j-1, k, 4, c]*vm1)
                   rhs[i, j, k, 5, c] += dy5ty1 *(
                         u[i, j+1, k, 5, c] - 2.0e0*u[i, j, k, 5, c] +
                          u[i, j-1, k, 5, c]) +
                         yycon3 * (qs[i, j+1, k, c] - 2.0e0*qs[i, j, k, c] +
                                   qs[i, j-1, k, c]) +
                         yycon4 * (vp1*vp1       - 2.0e0*vijk*vijk +
                                   vm1*vm1) +
                         yycon5 * (u[i, j+1, k, 5, c]*rho_i[i, j+1, k, c] -
                                   2.0e0*u[i, j, k, 5, c]*rho_i[i, j, k, c] +
                                   u[i, j-1, k, 5, c]*rho_i[i, j-1, k, c]) -
                         ty2 * ((c1*u[i, j+1, k, 5, c] -
                                 c2*square[i, j+1, k, c]) * vp1 -(
                                c1*u[i, j-1, k, 5, c] -
                                 c2*square[i, j-1, k, c]) * vm1)
                end
             end
          end

#---------------------------------------------------------------------
#         add fourth order eta-direction dissipation         
#---------------------------------------------------------------------
          if cell_start[2, c] > 0
             j = 1
             for m = 1:5
                for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
                   for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                      rhs[i, j, k, m, c] -= dssp *(5.0e0*u[i, j, k, m, c] - 4.0e0*u[i, j+1, k, m, c] + u[i, j+2, k, m, c])
                   end
                end
             end

             j = 2
             for m = 1:5
                for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
                   for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                      rhs[i, j, k, m, c] -= dssp *(-4.0e0*u[i, j-1, k, m, c] + 6.0e0*u[i, j, k, m, c] - 4.0e0*u[i, j+1, k, m, c] + u[i, j+2, k, m, c])
                   end
                end
             end
          end

          for m = 1:5
             for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
                for j = 3*cell_start[2, c]:cell_size[2, c]-3*cell_end[2, c]-1
                   for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                      rhs[i, j, k, m, c] -= dssp *(u[i, j-2, k, m, c] - 4.0e0*u[i, j-1, k, m, c] + 6.0*u[i, j, k, m, c] - 4.0e0*u[i, j+1, k, m, c] + u[i, j+2, k, m, c] )
                      #if j==0
                        #@info "j=-2 !!! rhs[$z][$i, $j, $k, $m, $c] =  $(rhs[i, j, k, m, c]) -- $(u[i, j-2, k, m, c])"
                      #end
                   end
                end
             end
          end

          if cell_end[2, c] > 0
             j = cell_size[2, c]-3
             for m = 1:5
                for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
                   for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                      rhs[i, j, k, m, c] -= dssp *(u[i, j-2, k, m, c] - 4.0e0*u[i, j-1, k, m, c] + 6.0e0*u[i, j, k, m, c] - 4.0e0*u[i, j+1, k, m, c] )
                   end
                end
             end

             j = cell_size[2, c]-2
             for m = 1:5
                for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
                   for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                      rhs[i, j, k, m, c] -= dssp *(u[i, j-2, k, m, c] - 4.0e0*u[i, j-1, k, m, c] + 5.0e0*u[i, j, k, m, c] )
                   end
                end
             end
          end


#---------------------------------------------------------------------
#         compute zeta-direction fluxes 
#---------------------------------------------------------------------
          for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
             for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
                for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                   wijk = ws[i, j, k, c]
                   wp1  = ws[i, j, k+1, c]
                   wm1  = ws[i, j, k-1, c]

                   rhs[i, j, k, 1, c] += dz1tz1 *(
                         u[i, j, k+1, 1, c] - 2.0e0*u[i, j, k, 1, c] +
                          u[i, j, k-1, 1, c]) -
                         tz2 * (u[i, j, k+1, 4, c] - u[i, j, k-1, 4, c])
                   rhs[i, j, k, 2, c] += dz2tz1 *(
                         u[i, j, k+1, 2, c] - 2.0e0*u[i, j, k, 2, c] +
                          u[i, j, k-1, 2, c]) +
                         zzcon2 * (us[i, j, k+1, c] - 2.0e0*us[i, j, k, c] +
                                   us[i, j, k-1, c]) -
                         tz2 * (u[i, j, k+1, 2, c]*wp1 -
                                u[i, j, k-1, 2, c]*wm1)
                   rhs[i, j, k, 3, c] += dz3tz1 *(
                         u[i, j, k+1, 3, c] - 2.0e0*u[i, j, k, 3, c] +
                          u[i, j, k-1, 3, c]) +
                         zzcon2 * (vs[i, j, k+1, c] - 2.0e0*vs[i, j, k, c] +
                                   vs[i, j, k-1, c]) -
                         tz2 * (u[i, j, k+1, 3, c]*wp1 -
                                u[i, j, k-1, 3, c]*wm1)
                   rhs[i, j, k, 4, c] += dz4tz1 *(
                         u[i, j, k+1, 4, c] - 2.0e0*u[i, j, k, 4, c] +
                          u[i, j, k-1, 4, c]) +
                         zzcon2*con43 * (wp1 - 2.0e0*wijk + wm1) -
                         tz2 * (u[i, j, k+1, 4, c]*wp1 -
                                u[i, j, k-1, 4, c]*wm1 +(
                                u[i, j, k+1, 5, c] - square[i, j, k+1, c] -
                                 u[i, j, k-1, 5, c] + square[i, j, k-1, c])*
                                c2)
                   rhs[i, j, k, 5, c] += dz5tz1 *(
                         u[i, j, k+1, 5, c] - 2.0e0*u[i, j, k, 5, c] +
                          u[i, j, k-1, 5, c]) +
                         zzcon3 * (qs[i, j, k+1, c] - 2.0e0*qs[i, j, k, c] + qs[i, j, k-1, c]) +
                         zzcon4 * (wp1*wp1 - 2.0e0*wijk*wijk + wm1*wm1) +
                         zzcon5 * (u[i, j, k+1, 5, c]*rho_i[i, j, k+1, c] -
                                   2.0e0*u[i, j, k, 5, c]*rho_i[i, j, k, c] +
                                   u[i, j, k-1, 5, c]*rho_i[i, j, k-1, c]) -
                         tz2 * ( (c1*u[i, j, k+1, 5, c] -
                                  c2*square[i, j, k+1, c])*wp1 -(
                                 c1*u[i, j, k-1, 5, c] -
                                  c2*square[i, j, k-1, c])*wm1)
                end
             end
          end

#---------------------------------------------------------------------
#         add fourth order zeta-direction dissipation                
#---------------------------------------------------------------------
          if cell_start[3, c] > 0
             k = 1
             for m = 1:5
                for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
                   for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                      rhs[i, j, k, m, c] -= dssp *(5.0e0*u[i, j, k, m, c] - 4.0e0*u[i, j, k+1, m, c] + u[i, j, k+2, m, c])
                   end
                end
             end

             k = 2
             for m = 1:5
                for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
                   for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                      rhs[i, j, k, m, c] -= dssp *(-4.0e0*u[i, j, k-1, m, c] + 6.0e0*u[i, j, k, m, c] - 4.0e0*u[i, j, k+1, m, c] + u[i, j, k+2, m, c])
                   end
                end
             end
          end

          for m = 1:5
             for k = 3*cell_start[3, c]:cell_size[3, c]-3*cell_end[3, c]-1
                for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
                   for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                      rhs[i, j, k, m, c]-= dssp *(u[i, j, k-2, m, c] - 4.0e0*u[i, j, k-1, m, c] + 6.0*u[i, j, k, m, c] - 4.0e0*u[i, j, k+1, m, c] + u[i, j, k+2, m, c])
                      #if k==0
                        #@info "k=-2 !!! rhs[$z][$i, $j, $k, $m, $c] =  $(rhs[i, j, k, m, c])"
                      #end
                   end
                end
             end
          end

          if cell_end[3, c] > 0
             k = cell_size[3, c]-3
             for m = 1:5
                for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
                   for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                      rhs[i, j, k, m, c] -= dssp *( u[i, j, k-2, m, c] - 4.0e0*u[i, j, k-1, m, c] + 6.0e0*u[i, j, k, m, c] - 4.0e0*u[i, j, k+1, m, c] )
                   end
                end
             end

             k = cell_size[3, c]-2
             for m = 1:5
                for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
                   for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                      rhs[i, j, k, m, c] -= dssp *(u[i, j, k-2, m, c] - 4.0e0*u[i, j, k-1, m, c] + 5.0e0*u[i, j, k, m, c] )
                   end
                end
             end
          end

          for m = 1:5
             for k = cell_start[3, c]:cell_size[3, c]-cell_end[3, c]-1
                for j = cell_start[2, c]:cell_size[2, c]-cell_end[2, c]-1
                   for i = cell_start[1, c]:cell_size[1, c]-cell_end[1, c]-1
                      rhs[i, j, k, m, c] *= dt
                   end
                end
             end
          end

       end


       if (timeron) timer_stop(t_rhs) end

       return nothing
end




