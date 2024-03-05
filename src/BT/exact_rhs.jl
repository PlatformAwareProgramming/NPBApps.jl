#---------------------------------------------------------------------
# compute the right hand side based on exact solution
#---------------------------------------------------------------------

function exact_rhs(z)

#---------------------------------------------------------------------
# loop over all cells owned by this node                   
#---------------------------------------------------------------------
       for c = 1:ncells

#---------------------------------------------------------------------
#         initialize                                  
#---------------------------------------------------------------------
          for m = 1:5
             for k = 0:cell_size[z][3, c]-1
                for j = 0:cell_size[z][2, c]-1
                   for i = 0:cell_size[z][1, c]-1
                      forcing[z][m, i, j, k, c] = 0.0e0
                   end
                end
             end
          end

#---------------------------------------------------------------------
# xi-direction flux differences                      
#---------------------------------------------------------------------
          for k = cell_start[z][3, c]:cell_size[z][3, c]-cell_end[z][3, c]-1
             zeta = float(k+cell_low[z][3, c]) * dnzm1
             for j = cell_start[z][2, c]:cell_size[z][2, c]-cell_end[z][2, c]-1
                eta = float(j+cell_low[z][2, c]) * dnym1

                for i = -2*(1-cell_start[z][1, c]):cell_size[z][1, c]+1-2*cell_end[z][1, c]
                   xi = float(i+cell_low[z][1, c]) * dnxm1

                   dtemp = exact_solution(xi, eta, zeta)
                   for m = 1:5
                      ue[z][i, m] = dtemp[m]
                   end

                   dtpp = 1.0e0 / dtemp[1]

                   for m = 2:5
                      buf[z][i, m] = dtpp * dtemp[m]
                   end

                   cuf[z][i]   = buf[z][i, 2] * buf[z][i, 2]
                   buf[z][i, 1] = cuf[z][i] + buf[z][i, 3] * buf[z][i, 3] + buf[z][i, 4] * buf[z][i, 4]
                   q[z][i] = 0.5e0*(buf[z][i, 2]*ue[z][i, 2] + buf[z][i, 3]*ue[z][i, 3] + buf[z][i, 4]*ue[z][i, 4])
                end

                for i = cell_start[z][1, c]:cell_size[z][1, c]-cell_end[z][1, c]-1
                   im1 = i-1
                   ip1 = i+1

                   forcing[z][1, i, j, k, c] = forcing[z][1, i, j, k, c] -
                       tx2*( ue[z][ip1, 2]-ue[z][im1, 2] )+
                       dx1tx1*(ue[z][ip1, 1]-2.0e0*ue[z][i, 1]+ue[z][im1, 1])

                   forcing[z][2, i, j, k, c] = forcing[z][2, i, j, k, c] - tx2 * ((
                      ue[z][ip1, 2]*buf[z][ip1, 2]+c2*(ue[z][ip1, 5]-q[z][ip1]))-(
                      ue[z][im1, 2]*buf[z][im1, 2]+c2*(ue[z][im1, 5]-q[z][im1])))+
                       xxcon1*(buf[z][ip1, 2]-2.0e0*buf[z][i, 2]+buf[z][im1, 2])+
                       dx2tx1*( ue[z][ip1, 2]-2.0e0* ue[z][i, 2]+ue[z][im1, 2])

                   forcing[z][3, i, j, k, c] = forcing[z][3, i, j, k, c] - tx2 * (
                       ue[z][ip1, 3]*buf[z][ip1, 2]-ue[z][im1, 3]*buf[z][im1, 2])+
                       xxcon2*(buf[z][ip1, 3]-2.0e0*buf[z][i, 3]+buf[z][im1, 3])+
                       dx3tx1*( ue[z][ip1, 3]-2.0e0*ue[z][i, 3] +ue[z][im1, 3])

                   forcing[z][4, i, j, k, c] = forcing[z][4, i, j, k, c] - tx2*(
                       ue[z][ip1, 4]*buf[z][ip1, 2]-ue[z][im1, 4]*buf[z][im1, 2])+
                       xxcon2*(buf[z][ip1, 4]-2.0e0*buf[z][i, 4]+buf[z][im1, 4])+
                       dx4tx1*( ue[z][ip1, 4]-2.0e0* ue[z][i, 4]+ ue[z][im1, 4])

                   forcing[z][5, i, j, k, c] = forcing[z][5, i, j, k, c] - tx2*(
                       buf[z][ip1, 2]*(c1*ue[z][ip1, 5]-c2*q[z][ip1])-
                       buf[z][im1, 2]*(c1*ue[z][im1, 5]-c2*q[z][im1]))+
                       0.5e0*xxcon3*(buf[z][ip1, 1]-2.0e0*buf[z][i, 1]+ buf[z][im1, 1])+
                       xxcon4*(cuf[z][ip1]-2.0e0*cuf[z][i]+cuf[z][im1])+
                       xxcon5*(buf[z][ip1, 5]-2.0e0*buf[z][i, 5]+buf[z][im1, 5])+
                       dx5tx1*( ue[z][ip1, 5]-2.0e0* ue[z][i, 5]+ ue[z][im1, 5])
                end

#---------------------------------------------------------------------
# Fourth-order dissipation                         
#---------------------------------------------------------------------
                if cell_start[z][1, c] > 0
                   for m = 1:5
                      i = 1
                      forcing[z][m, i, j, k, c] = forcing[z][m, i, j, k, c] - dssp *(
                          5.0e0*ue[z][i, m] - 4.0e0*ue[z][i+1, m] +ue[z][i+2, m])
                      i = 2
                      forcing[z][m, i, j, k, c] = forcing[z][m, i, j, k, c] - dssp *(
                         -4.0e0*ue[z][i-1, m] + 6.0e0*ue[z][i, m] -
                           4.0e0*ue[z][i+1, m] +       ue[z][i+2, m])
                   end
                end

                for m = 1:5
                   for i = cell_start[z][1, c]*3:cell_size[z][1, c]-3*cell_end[z][1, c]-1
                      forcing[z][m, i, j, k, c] = forcing[z][m, i, j, k, c] - dssp*(
                         ue[z][i-2, m] - 4.0e0*ue[z][i-1, m] +
                          6.0e0*ue[z][i, m] - 4.0e0*ue[z][i+1, m] + ue[z][i+2, m])
                   end
                end

                if cell_end[z][1, c] > 0
                   for m = 1:5
                      i = cell_size[z][1, c]-3
                      forcing[z][m, i, j, k, c] = forcing[z][m, i, j, k, c] - dssp *(
                         ue[z][i-2, m] - 4.0e0*ue[z][i-1, m] +
                          6.0e0*ue[z][i, m] - 4.0e0*ue[z][i+1, m])
                      i = cell_size[z][1, c]-2
                      forcing[z][m, i, j, k, c] = forcing[z][m, i, j, k, c] - dssp *(
                         ue[z][i-2, m] - 4.0e0*ue[z][i-1, m] + 5.0e0*ue[z][i, m])
                   end
                end

             end
          end
#---------------------------------------------------------------------
#  eta-direction flux differences             
#---------------------------------------------------------------------
          for k = cell_start[z][3, c]:cell_size[z][3, c]-cell_end[z][3, c]-1
             zeta = float(k+cell_low[z][3, c]) * dnzm1
             for i = cell_start[z][1, c]:cell_size[z][1, c]-cell_end[z][1, c]-1
                xi = float(i+cell_low[z][1, c]) * dnxm1

                for j = -2*(1-cell_start[z][2, c]):cell_size[z][2, c]+1-2*cell_end[z][2, c]
                   eta = float(j+cell_low[z][2, c]) * dnym1

                   dtemp = exact_solution(xi, eta, zeta)
                   for m = 1:5
                      ue[z][j, m] = dtemp[m]
                   end
                   dtpp = 1.0e0/dtemp[1]

                   for m = 2:5
                      buf[z][j, m] = dtpp * dtemp[m]
                   end

                   cuf[z][j]   = buf[z][j, 3] * buf[z][j, 3]
                   buf[z][j, 1] = cuf[z][j] + buf[z][j, 2] * buf[z][j, 2] + buf[z][j, 4] * buf[z][j, 4]
                   q[z][j] = 0.5e0*(buf[z][j, 2]*ue[z][j, 2] + buf[z][j, 3]*ue[z][j, 3] + buf[z][j, 4]*ue[z][j, 4])
                end

                for j = cell_start[z][2, c]:cell_size[z][2, c]-cell_end[z][2, c]-1
                   jm1 = j-1
                   jp1 = j+1

                   forcing[z][1, i, j, k, c] = forcing[z][1, i, j, k, c] -
                      ty2*( ue[z][jp1, 3]-ue[z][jm1, 3] )+
                      dy1ty1*(ue[z][jp1, 1]-2.0e0*ue[z][j, 1]+ue[z][jm1, 1])

                   forcing[z][2, i, j, k, c] = forcing[z][2, i, j, k, c] - ty2*(
                      ue[z][jp1, 2]*buf[z][jp1, 3]-ue[z][jm1, 2]*buf[z][jm1, 3])+
                      yycon2*(buf[z][jp1, 2]-2.0e0*buf[z][j, 2]+buf[z][jm1, 2])+
                      dy2ty1*( ue[z][jp1, 2]-2.0* ue[z][j, 2]+ ue[z][jm1, 2])

                   forcing[z][3, i, j, k, c] = forcing[z][3, i, j, k, c] - ty2*((
                      ue[z][jp1, 3]*buf[z][jp1, 3]+c2*(ue[z][jp1, 5]-q[z][jp1]))-(
                      ue[z][jm1, 3]*buf[z][jm1, 3]+c2*(ue[z][jm1, 5]-q[z][jm1])))+
                      yycon1*(buf[z][jp1, 3]-2.0e0*buf[z][j, 3]+buf[z][jm1, 3])+
                      dy3ty1*( ue[z][jp1, 3]-2.0e0*ue[z][j, 3] +ue[z][jm1, 3])

                   forcing[z][4, i, j, k, c] = forcing[z][4, i, j, k, c] - ty2*(
                      ue[z][jp1, 4]*buf[z][jp1, 3]-ue[z][jm1, 4]*buf[z][jm1, 3])+
                      yycon2*(buf[z][jp1, 4]-2.0e0*buf[z][j, 4]+buf[z][jm1, 4])+
                      dy4ty1*( ue[z][jp1, 4]-2.0e0*ue[z][j, 4]+ ue[z][jm1, 4])

                   forcing[z][5, i, j, k, c] = forcing[z][5, i, j, k, c] - ty2*(
                      buf[z][jp1, 3]*(c1*ue[z][jp1, 5]-c2*q[z][jp1])-
                      buf[z][jm1, 3]*(c1*ue[z][jm1, 5]-c2*q[z][jm1]))+
                      0.5e0*yycon3*(buf[z][jp1, 1]-2.0e0*buf[z][j, 1]+
                                    buf[z][jm1, 1])+
                      yycon4*(cuf[z][jp1]-2.0e0*cuf[z][j]+cuf[z][jm1])+
                      yycon5*(buf[z][jp1, 5]-2.0e0*buf[z][j, 5]+buf[z][jm1, 5])+
                      dy5ty1*(ue[z][jp1, 5]-2.0e0*ue[z][j, 5]+ue[z][jm1, 5])
                end

#---------------------------------------------------------------------
# Fourth-order dissipation                      
#---------------------------------------------------------------------
                if cell_start[z][2, c] > 0
                   for m = 1:5
                      j = 1
                      forcing[z][m, i, j, k, c] = forcing[z][m, i, j, k, c] - dssp *(
                          5.0e0*ue[z][j, m] - 4.0e0*ue[z][j+1, m] +ue[z][j+2, m])
                      j = 2
                      forcing[z][m, i, j, k, c] = forcing[z][m, i, j, k, c] - dssp *(
                         -4.0e0*ue[z][j-1, m] + 6.0e0*ue[z][j, m] -
                           4.0e0*ue[z][j+1, m] +       ue[z][j+2, m])
                   end
                end

                for m = 1:5
                   for j = cell_start[z][2, c]*3:cell_size[z][2, c]-3*cell_end[z][2, c]-1
                      forcing[z][m, i, j, k, c] = forcing[z][m, i, j, k, c] - dssp*(
                         ue[z][j-2, m] - 4.0e0*ue[z][j-1, m] +
                          6.0e0*ue[z][j, m] - 4.0e0*ue[z][j+1, m] + ue[z][j+2, m])
                   end
                end
                if cell_end[z][2, c] > 0
                   for m = 1:5
                      j = cell_size[z][2, c]-3
                      forcing[z][m, i, j, k, c] = forcing[z][m, i, j, k, c] - dssp *(
                         ue[z][j-2, m] - 4.0e0*ue[z][j-1, m] +
                          6.0e0*ue[z][j, m] - 4.0e0*ue[z][j+1, m])
                      j = cell_size[z][2, c]-2
                      forcing[z][m, i, j, k, c] = forcing[z][m, i, j, k, c] - dssp *(
                         ue[z][j-2, m] - 4.0e0*ue[z][j-1, m] + 5.0e0*ue[z][j, m])

                   end
                end

             end
          end

#---------------------------------------------------------------------
# zeta-direction flux differences                      
#---------------------------------------------------------------------
          for j = cell_start[z][2, c]:cell_size[z][2, c]-cell_end[z][2, c]-1
             eta = float(j+cell_low[z][2, c]) * dnym1
             for i = cell_start[z][1, c]:cell_size[z][1, c]-cell_end[z][1, c]-1
                xi = float(i+cell_low[z][1, c]) * dnxm1

                for k = -2*(1-cell_start[z][3, c]):cell_size[z][3, c]+1-2*cell_end[z][3, c]
                   zeta = float(k+cell_low[z][3, c]) * dnzm1

                   dtemp = exact_solution(xi, eta, zeta)
                   for m = 1:5
                      ue[z][k, m] = dtemp[m]
                   end

                   dtpp = 1.0e0/dtemp[1]

                   for m = 2:5
                      buf[z][k, m] = dtpp * dtemp[m]
                   end

                   cuf[z][k]   = buf[z][k, 4] * buf[z][k, 4]
                   buf[z][k, 1] = cuf[z][k] + buf[z][k, 2] * buf[z][k, 2] +
                              buf[z][k, 3] * buf[z][k, 3]
                   q[z][k] = 0.5e0*(buf[z][k, 2]*ue[z][k, 2] + buf[z][k, 3]*ue[z][k, 3] +
                                 buf[z][k, 4]*ue[z][k, 4])
                end

                for k = cell_start[z][3, c]:cell_size[z][3, c]-cell_end[z][3, c]-1
                   km1 = k-1
                   kp1 = k+1

                   forcing[z][1, i, j, k, c] = forcing[z][1, i, j, k, c] -
                       tz2*( ue[z][kp1, 4]-ue[z][km1, 4] )+
                       dz1tz1*(ue[z][kp1, 1]-2.0e0*ue[z][k, 1]+ue[z][km1, 1])

                   forcing[z][2, i, j, k, c] = forcing[z][2, i, j, k, c] - tz2 * (
                       ue[z][kp1, 2]*buf[z][kp1, 4]-ue[z][km1, 2]*buf[z][km1, 4])+
                       zzcon2*(buf[z][kp1, 2]-2.0e0*buf[z][k, 2]+buf[z][km1, 2])+
                       dz2tz1*( ue[z][kp1, 2]-2.0e0* ue[z][k, 2]+ ue[z][km1, 2])

                   forcing[z][3, i, j, k, c] = forcing[z][3, i, j, k, c] - tz2 * (
                       ue[z][kp1, 3]*buf[z][kp1, 4]-ue[z][km1, 3]*buf[z][km1, 4])+
                       zzcon2*(buf[z][kp1, 3]-2.0e0*buf[z][k, 3]+buf[z][km1, 3])+
                       dz3tz1*(ue[z][kp1, 3]-2.0e0*ue[z][k, 3]+ue[z][km1, 3])

                   forcing[z][4, i, j, k, c] = forcing[z][4, i, j, k, c] - tz2 * ((
                      ue[z][kp1, 4]*buf[z][kp1, 4]+c2*(ue[z][kp1, 5]-q[z][kp1]))-(
                      ue[z][km1, 4]*buf[z][km1, 4]+c2*(ue[z][km1, 5]-q[z][km1])))+
                      zzcon1*(buf[z][kp1, 4]-2.0e0*buf[z][k, 4]+buf[z][km1, 4])+
                      dz4tz1*( ue[z][kp1, 4]-2.0e0*ue[z][k, 4] +ue[z][km1, 4])

                   forcing[z][5, i, j, k, c] = forcing[z][5, i, j, k, c] - tz2 * (
                       buf[z][kp1, 4]*(c1*ue[z][kp1, 5]-c2*q[z][kp1])-
                       buf[z][km1, 4]*(c1*ue[z][km1, 5]-c2*q[z][km1]))+
                       0.5e0*zzcon3*(buf[z][kp1, 1]-2.0e0*buf[z][k, 1]+
                                    buf[z][km1, 1])+
                       zzcon4*(cuf[z][kp1]-2.0e0*cuf[z][k]+cuf[z][km1])+
                       zzcon5*(buf[z][kp1, 5]-2.0e0*buf[z][k, 5]+buf[z][km1, 5])+
                       dz5tz1*( ue[z][kp1, 5]-2.0e0*ue[z][k, 5]+ ue[z][km1, 5])
                end

#---------------------------------------------------------------------
# Fourth-order dissipation                        
#---------------------------------------------------------------------
                if cell_start[z][3, c] > 0
                   for m = 1:5
                      k = 1
                      forcing[z][m, i, j, k, c] = forcing[z][m, i, j, k, c] - dssp *(
                          5.0e0*ue[z][k, m] - 4.0e0*ue[z][k+1, m] +ue[z][k+2, m])
                      k = 2
                      forcing[z][m, i, j, k, c] = forcing[z][m, i, j, k, c] - dssp *(
                         -4.0e0*ue[z][k-1, m] + 6.0e0*ue[z][k, m] -
                           4.0e0*ue[z][k+1, m] +       ue[z][k+2, m])
                   end
                end

                for m = 1:5
                   for k = cell_start[z][3, c]*3:cell_size[z][3, c]-3*cell_end[z][3, c]-1
                      forcing[z][m, i, j, k, c] = forcing[z][m, i, j, k, c] - dssp*(
                         ue[z][k-2, m] - 4.0e0*ue[z][k-1, m] +
                          6.0e0*ue[z][k, m] - 4.0e0*ue[z][k+1, m] + ue[z][k+2, m])
                   end
                end

                if cell_end[z][3, c] > 0
                   for m = 1:5
                      k = cell_size[z][3, c]-3
                      forcing[z][m, i, j, k, c] = forcing[z][m, i, j, k, c] - dssp *(
                         ue[z][k-2, m] - 4.0e0*ue[z][k-1, m] +
                          6.0e0*ue[z][k, m] - 4.0e0*ue[z][k+1, m])
                      k = cell_size[z][3, c]-2
                      forcing[z][m, i, j, k, c] = forcing[z][m, i, j, k, c] - dssp *(
                         ue[z][k-2, m] - 4.0e0*ue[z][k-1, m] + 5.0e0*ue[z][k, m])
                   end
                end

             end
          end
#---------------------------------------------------------------------
# now change the sign of the forcing function, 
#---------------------------------------------------------------------
          for m = 1:5
             for k = cell_start[z][3, c]:cell_size[z][3, c]-cell_end[z][3, c]-1
                for j = cell_start[z][2, c]:cell_size[z][2, c]-cell_end[z][2, c]-1
                   for i = cell_start[z][1, c]:cell_size[z][1, c]-cell_end[z][1, c]-1
                      forcing[z][m, i, j, k, c] = -1.0e0 * forcing[z][m, i, j, k, c]
                   end
                end
             end
          end

#---------------------------------------------------------------------
#      cell loop
#---------------------------------------------------------------------
       end

       return nothing
end





