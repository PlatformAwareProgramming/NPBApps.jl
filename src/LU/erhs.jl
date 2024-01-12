#---------------------------------------------------------------------
#
#   compute the right hand side based on exact solution
#
#---------------------------------------------------------------------

 function erhs(nx0, ny0, nz0)

      dsspm = dssp

      for k = 1:nz
         for j = 1:ny
            for i = 1:nx
               for m = 1:5
                  frct[ m, i, j, k ] = 0.0e+00
               end
            end
         end
      end

      for k = 1:nz
         zeta = ( float(k-1) ) / ( nz - 1 )
         for j = 1:ny
            jglob = jpt + j
            eta = ( float(jglob-1) ) / ( ny0 - 1 )
            for i = 1:nx
               iglob = ipt + i
               xi = ( float(iglob-1) ) / ( nx0 - 1 )
               for m = 1:5
                  rsd[m, i, j, k] =  ce[m, 1]+
                        ce[m, 2] * xi+
                        ce[m, 3] * eta+
                        ce[m, 4] * zeta+
                        ce[m, 5] * xi * xi+
                        ce[m, 6] * eta * eta+
                        ce[m, 7] * zeta * zeta+
                        ce[m, 8] * xi * xi * xi+
                        ce[m, 9] * eta * eta * eta+
                        ce[m, 10] * zeta * zeta * zeta+
                        ce[m, 11] * xi * xi * xi * xi+
                        ce[m, 12] * eta * eta * eta * eta+
                        ce[m, 13] * zeta * zeta * zeta * zeta
               end
            end
         end
      end

#---------------------------------------------------------------------
#   xi-direction flux differences
#---------------------------------------------------------------------
#
#   iex = flag : iex = 0  north/south communication
#              : iex = 1  east/west communication
#
#---------------------------------------------------------------------
      iex   = 0

#---------------------------------------------------------------------
#   communicate and receive/send two rows of data
#---------------------------------------------------------------------
      exchange_3(rsd, iex,
                  comm_solve, 
                  buf,
                  buf1,
                  south,
                  east,
                  north,
                  west,
                  nx,
                  ny,
                  nz,
               )

      L1 = 0
      if (north == -1) L1 = 1 end
      L2 = nx + 1
      if (south == -1) L2 = nx end

      ist1 = 1
      iend1 = nx
      if (north == -1) ist1 = 4 end
      if (south == -1) iend1 = nx - 3 end

      for k = 2:nz - 1
         for j = jst:jend
            for i = L1:L2
               flux[1, i, j, k] = rsd[2, i, j, k]
               u21 = rsd[2, i, j, k] / rsd[1, i, j, k]
               q = 0.50e+00 * (  rsd[2, i, j, k] * rsd[2, i, j, k]+
                                rsd[3, i, j, k] * rsd[3, i, j, k]+
                                rsd[4, i, j, k] * rsd[4, i, j, k] )/
                             rsd[1, i, j, k]
               flux[2, i, j, k] = rsd[2, i, j, k] * u21 + c2 *(
                                rsd[5, i, j, k] - q )
               flux[3, i, j, k] = rsd[3, i, j, k] * u21
               flux[4, i, j, k] = rsd[4, i, j, k] * u21
               flux[5, i, j, k] = ( c1 * rsd[5, i, j, k] - c2 * q ) * u21
            end
         end
      end

      for k = 2:nz - 1
         for j = jst:jend
            for i = ist:iend
               for m = 1:5
                  frct[m, i, j, k] =  frct[m, i, j, k]-
                          tx2 * ( flux[m, i+1, j, k] - flux[m, i-1, j, k] )
               end
            end
            for i = ist:L2
               tmp = 1.0e+00 / rsd[1, i, j, k]

               u21i = tmp * rsd[2, i, j, k]
               u31i = tmp * rsd[3, i, j, k]
               u41i = tmp * rsd[4, i, j, k]
               u51i = tmp * rsd[5, i, j, k]

               tmp = 1.0e+00 / rsd[1, i-1, j, k]

               u21im1 = tmp * rsd[2, i-1, j, k]
               u31im1 = tmp * rsd[3, i-1, j, k]
               u41im1 = tmp * rsd[4, i-1, j, k]
               u51im1 = tmp * rsd[5, i-1, j, k]

               flux[2, i, j, k] = (4.0e+00/3.0e+00) * tx3 *(
                               u21i - u21im1 )
               flux[3, i, j, k] = tx3 * ( u31i - u31im1 )
               flux[4, i, j, k] = tx3 * ( u41i - u41im1 )
               flux[5, i, j, k] = 0.50e+00 * ( 1.0e+00 - c1*c5 )*
                     tx3 * ( ( u21i  ^2 + u31i  ^2 + u41i  ^2 )-
                             ( u21im1^2 + u31im1^2 + u41im1^2 ) )+
                     (1.0e+00/6.0e+00)*
                     tx3 * ( u21i^2 - u21im1^2 )+
                     c1 * c5 * tx3 * ( u51i - u51im1 )
            end

            for i = ist:iend
               frct[1, i, j, k] = frct[1, i, j, k]+
                     dx1 * tx1 * (            rsd[1, i-1, j, k]-
                                    2.0e+00 * rsd[1, i, j, k]+
                                              rsd[1, i+1, j, k] )
               frct[2, i, j, k] = frct[2, i, j, k]+
                  tx3 * c3 * c4 * ( flux[2, i+1, j, k] - flux[2, i, j, k] )+
                     dx2 * tx1 * (            rsd[2, i-1, j, k]-
                                    2.0e+00 * rsd[2, i, j, k]+
                                              rsd[2, i+1, j, k] )
               frct[3, i, j, k] = frct[3, i, j, k]+
                  tx3 * c3 * c4 * ( flux[3, i+1, j, k] - flux[3, i, j, k] )+
                     dx3 * tx1 * (            rsd[3, i-1, j, k]-
                                    2.0e+00 * rsd[3, i, j, k]+
                                              rsd[3, i+1, j, k] )
               frct[4, i, j, k] = frct[4, i, j, k]+
                   tx3 * c3 * c4 * ( flux[4, i+1, j, k] - flux[4, i, j, k] )+
                     dx4 * tx1 * (            rsd[4, i-1, j, k]-
                                    2.0e+00 * rsd[4, i, j, k]+
                                              rsd[4, i+1, j, k] )
               frct[5, i, j, k] = frct[5, i, j, k]+
                  tx3 * c3 * c4 * ( flux[5, i+1, j, k] - flux[5, i, j, k] )+
                     dx5 * tx1 * (            rsd[5, i-1, j, k]-
                                    2.0e+00 * rsd[5, i, j, k]+
                                              rsd[5, i+1, j, k] )
            end

#---------------------------------------------------------------------
#   Fourth-order dissipation
#---------------------------------------------------------------------
            if north == -1
             for m = 1:5
               frct[m, 2, j, k] = frct[m, 2, j, k]-
                  dsspm * ( + 5.0e+00 * rsd[m, 2, j, k]-
                              4.0e+00 * rsd[m, 3, j, k]+
                                        rsd[m, 4, j, k] )
               frct[m, 3, j, k] = frct[m, 3, j, k]-
                  dsspm * ( - 4.0e+00 * rsd[m, 2, j, k]+
                              6.0e+00 * rsd[m, 3, j, k]-
                              4.0e+00 * rsd[m, 4, j, k]+
                                        rsd[m, 5, j, k] )
             end
            end

            for i = ist1:iend1
               for m = 1:5
                  frct[m, i, j, k] = frct[m, i, j, k]-
                     dsspm * (            rsd[m, i-2, j, k]-
                                4.0e+00 * rsd[m, i-1, j, k]+
                                6.0e+00 * rsd[m, i, j, k]-
                                4.0e+00 * rsd[m, i+1, j, k]+
                                          rsd[m, i+2, j, k] )
               end
            end

            if south == -1
             for m = 1:5
               frct[m, nx-2, j, k] = frct[m, nx-2, j, k]-
                  dsspm * (             rsd[m, nx-4, j, k]-
                              4.0e+00 * rsd[m, nx-3, j, k]+
                              6.0e+00 * rsd[m, nx-2, j, k]-
                              4.0e+00 * rsd[m, nx-1, j, k]  )
               frct[m, nx-1, j, k] = frct[m, nx-1, j, k]-
                  dsspm * (             rsd[m, nx-3, j, k]-
                              4.0e+00 * rsd[m, nx-2, j, k]+
                              5.0e+00 * rsd[m, nx-1, j, k] )
             end
            end

         end
      end

#---------------------------------------------------------------------
#   eta-direction flux differences
#---------------------------------------------------------------------
#
#   iex = flag : iex = 0  north/south communication
#              : iex = 1  east/west communication
#
#---------------------------------------------------------------------
      iex   = 1

#---------------------------------------------------------------------
#   communicate and receive/send two rows of data
#---------------------------------------------------------------------
      exchange_3(rsd, iex,
                  comm_solve, 
                  buf,
                  buf1,
                  south,
                  east,
                  north,
                  west,
                  nx,
                  ny,
                  nz,
               )

      L1 = 0
      if (west == -1) L1 = 1 end
      L2 = ny + 1
      if (east == -1) L2 = ny end

      jst1 = 1
      jend1 = ny
      if (west == -1) jst1 = 4 end
      if (east == -1) jend1 = ny - 3 end

      for k = 2:nz - 1
         for j = L1:L2
            for i = ist:iend
               flux[1, i, j, k] = rsd[3, i, j, k]
               u31 = rsd[3, i, j, k] / rsd[1, i, j, k]
               q = 0.50e+00 * (  rsd[2, i, j, k] * rsd[2, i, j, k]+
                                rsd[3, i, j, k] * rsd[3, i, j, k]+
                                rsd[4, i, j, k] * rsd[4, i, j, k] )/
                             rsd[1, i, j, k]
               flux[2, i, j, k] = rsd[2, i, j, k] * u31
               flux[3, i, j, k] = rsd[3, i, j, k] * u31 + c2 *(
                              rsd[5, i, j, k] - q )
               flux[4, i, j, k] = rsd[4, i, j, k] * u31
               flux[5, i, j, k] = ( c1 * rsd[5, i, j, k] - c2 * q ) * u31
            end
         end
      end

      for k = 2:nz - 1
         for i = ist:iend
            for j = jst:jend
               for m = 1:5
                  frct[m, i, j, k] =  frct[m, i, j, k]-
                        ty2 * ( flux[m, i, j+1, k] - flux[m, i, j-1, k] )
               end
            end
         end

         for j = jst:L2
            for i = ist:iend
               tmp = 1.0e+00 / rsd[1, i, j, k]

               u21j = tmp * rsd[2, i, j, k]
               u31j = tmp * rsd[3, i, j, k]
               u41j = tmp * rsd[4, i, j, k]
               u51j = tmp * rsd[5, i, j, k]

               tmp = 1.0e+00 / rsd[1, i, j-1, k]

               u21jm1 = tmp * rsd[2, i, j-1, k]
               u31jm1 = tmp * rsd[3, i, j-1, k]
               u41jm1 = tmp * rsd[4, i, j-1, k]
               u51jm1 = tmp * rsd[5, i, j-1, k]

               flux[2, i, j, k] = ty3 * ( u21j - u21jm1 )
               flux[3, i, j, k] = (4.0e+00/3.0e+00) * ty3 *(
                              u31j - u31jm1 )
               flux[4, i, j, k] = ty3 * ( u41j - u41jm1 )
               flux[5, i, j, k] = 0.50e+00 * ( 1.0e+00 - c1*c5 )*
                     ty3 * ( ( u21j  ^2 + u31j  ^2 + u41j  ^2 )-
                             ( u21jm1^2 + u31jm1^2 + u41jm1^2 ) )+
                     (1.0e+00/6.0e+00)*
                     ty3 * ( u31j^2 - u31jm1^2 )+
                     c1 * c5 * ty3 * ( u51j - u51jm1 )
            end
         end

         for j = jst:jend
            for i = ist:iend
               frct[1, i, j, k] = frct[1, i, j, k]+
                     dy1 * ty1 * (            rsd[1, i, j-1, k]-
                                    2.0e+00 * rsd[1, i, j, k]+
                                              rsd[1, i, j+1, k] )
               frct[2, i, j, k] = frct[2, i, j, k]+
                 ty3 * c3 * c4 * ( flux[2, i, j+1, k] - flux[2, i, j, k] )+
                     dy2 * ty1 * (            rsd[2, i, j-1, k]-
                                    2.0e+00 * rsd[2, i, j, k]+
                                              rsd[2, i, j+1, k] )
               frct[3, i, j, k] = frct[3, i, j, k]+
                 ty3 * c3 * c4 * ( flux[3, i, j+1, k] - flux[3, i, j, k] )+
                     dy3 * ty1 * (            rsd[3, i, j-1, k]-
                                    2.0e+00 * rsd[3, i, j, k]+
                                              rsd[3, i, j+1, k] )
               frct[4, i, j, k] = frct[4, i, j, k]+
                 ty3 * c3 * c4 * ( flux[4, i, j+1, k] - flux[4, i, j, k] )+
                     dy4 * ty1 * (            rsd[4, i, j-1, k]-
                                    2.0e+00 * rsd[4, i, j, k]+
                                              rsd[4, i, j+1, k] )
               frct[5, i, j, k] = frct[5, i, j, k]+
                 ty3 * c3 * c4 * ( flux[5, i, j+1, k] - flux[5, i, j, k] )+
                     dy5 * ty1 * (            rsd[5, i, j-1, k]-
                                    2.0e+00 * rsd[5, i, j, k]+
                                              rsd[5, i, j+1, k] )
            end
         end

#---------------------------------------------------------------------
#   fourth-order dissipation
#---------------------------------------------------------------------
         if west == -1
            for i = ist:iend
             for m = 1:5
               frct[m, i, 2, k] = frct[m, i, 2, k]-
                  dsspm * ( + 5.0e+00 * rsd[m, i, 2, k]-
                              4.0e+00 * rsd[m, i, 3, k]+
                                        rsd[m, i, 4, k] )
               frct[m, i, 3, k] = frct[m, i, 3, k]-
                  dsspm * ( - 4.0e+00 * rsd[m, i, 2, k]+
                              6.0e+00 * rsd[m, i, 3, k]-
                              4.0e+00 * rsd[m, i, 4, k]+
                                        rsd[m, i, 5, k] )
             end
            end
         end

         for j = jst1:jend1
            for i = ist:iend
               for m = 1:5
                  frct[m, i, j, k] = frct[m, i, j, k]-
                     dsspm * (            rsd[m, i, j-2, k]-
                               4.0e+00 * rsd[m, i, j-1, k]+
                               6.0e+00 * rsd[m, i, j, k]-
                               4.0e+00 * rsd[m, i, j+1, k]+
                                         rsd[m, i, j+2, k] )
               end
            end
         end

         if east == -1
            for i = ist:iend
             for m = 1:5
               frct[m, i, ny-2, k] = frct[m, i, ny-2, k]-
                  dsspm * (             rsd[m, i, ny-4, k]-
                              4.0e+00 * rsd[m, i, ny-3, k]+
                              6.0e+00 * rsd[m, i, ny-2, k]-
                              4.0e+00 * rsd[m, i, ny-1, k]  )
               frct[m, i, ny-1, k] = frct[m, i, ny-1, k]-
                  dsspm * (             rsd[m, i, ny-3, k]-
                              4.0e+00 * rsd[m, i, ny-2, k]+
                              5.0e+00 * rsd[m, i, ny-1, k]  )
             end
            end
         end

      end

#---------------------------------------------------------------------
#   zeta-direction flux differences
#---------------------------------------------------------------------
      for k = 1:nz
         for j = jst:jend
            for i = ist:iend
               flux[1, i, j, k] = rsd[4, i, j, k]
               u41 = rsd[4, i, j, k] / rsd[1, i, j, k]
               q = 0.50e+00 * (  rsd[2, i, j, k] * rsd[2, i, j, k]+
                                rsd[3, i, j, k] * rsd[3, i, j, k]+
                                rsd[4, i, j, k] * rsd[4, i, j, k] )/
                             rsd[1, i, j, k]
               flux[2, i, j, k] = rsd[2, i, j, k] * u41
               flux[3, i, j, k] = rsd[3, i, j, k] * u41
               flux[4, i, j, k] = rsd[4, i, j, k] * u41 + c2 *(
                                rsd[5, i, j, k] - q )
               flux[5, i, j, k] = ( c1 * rsd[5, i, j, k] - c2 * q ) * u41
            end
         end
      end

      for k = 2:nz - 1
         for j = jst:jend
            for i = ist:iend
               for m = 1:5
                  frct[m, i, j, k] =  frct[m, i, j, k]-
                         tz2 * ( flux[m, i, j, k+1] - flux[m, i, j, k-1] )
               end
            end
         end
      end

      for k = 2:nz
         for j = jst:jend
            for i = ist:iend
               tmp = 1.0e+00 / rsd[1, i, j, k]

               u21k = tmp * rsd[2, i, j, k]
               u31k = tmp * rsd[3, i, j, k]
               u41k = tmp * rsd[4, i, j, k]
               u51k = tmp * rsd[5, i, j, k]

               tmp = 1.0e+00 / rsd[1, i, j, k-1]

               u21km1 = tmp * rsd[2, i, j, k-1]
               u31km1 = tmp * rsd[3, i, j, k-1]
               u41km1 = tmp * rsd[4, i, j, k-1]
               u51km1 = tmp * rsd[5, i, j, k-1]

               flux[2, i, j, k] = tz3 * ( u21k - u21km1 )
               flux[3, i, j, k] = tz3 * ( u31k - u31km1 )
               flux[4, i, j, k] = (4.0e+00/3.0e+00) * tz3 * ( u41k-
                              u41km1 )
               flux[5, i, j, k] = 0.50e+00 * ( 1.0e+00 - c1*c5 )*
                     tz3 * ( ( u21k  ^2 + u31k  ^2 + u41k  ^2 )-
                             ( u21km1^2 + u31km1^2 + u41km1^2 ) )+
                     (1.0e+00/6.0e+00)*
                     tz3 * ( u41k^2 - u41km1^2 )+
                     c1 * c5 * tz3 * ( u51k - u51km1 )
            end
         end
      end

      for k = 2:nz - 1
         for j = jst:jend
            for i = ist:iend
               frct[1, i, j, k] = frct[1, i, j, k]+
                     dz1 * tz1 * (            rsd[1, i, j, k+1]-
                                    2.0e+00 * rsd[1, i, j, k]+
                                              rsd[1, i, j, k-1] )
               frct[2, i, j, k] = frct[2, i, j, k]+
                 tz3 * c3 * c4 * ( flux[2, i, j, k+1] - flux[2, i, j, k] )+
                     dz2 * tz1 * (            rsd[2, i, j, k+1]-
                                    2.0e+00 * rsd[2, i, j, k]+
                                              rsd[2, i, j, k-1] )
               frct[3, i, j, k] = frct[3, i, j, k]+
                 tz3 * c3 * c4 * ( flux[3, i, j, k+1] - flux[3, i, j, k] )+
                     dz3 * tz1 * (            rsd[3, i, j, k+1]-
                                    2.0e+00 * rsd[3, i, j, k]+
                                              rsd[3, i, j, k-1] )
               frct[4, i, j, k] = frct[4, i, j, k]+
                 tz3 * c3 * c4 * ( flux[4, i, j, k+1] - flux[4, i, j, k] )+
                     dz4 * tz1 * (            rsd[4, i, j, k+1]-
                                    2.0e+00 * rsd[4, i, j, k]+
                                              rsd[4, i, j, k-1] )
               frct[5, i, j, k] = frct[5, i, j, k]+
                 tz3 * c3 * c4 * ( flux[5, i, j, k+1] - flux[5, i, j, k] )+
                     dz5 * tz1 * (            rsd[5, i, j, k+1]-
                                    2.0e+00 * rsd[5, i, j, k]+
                                              rsd[5, i, j, k-1] )
            end
         end
      end

#---------------------------------------------------------------------
#   fourth-order dissipation
#---------------------------------------------------------------------
      for j = jst:jend
         for i = ist:iend
            for m = 1:5
               frct[m, i, j, 2] = frct[m, i, j, 2]-
                  dsspm * ( + 5.0e+00 * rsd[m, i, j, 2]-
                              4.0e+00 * rsd[m, i, j, 3]+
                                        rsd[m, i, j, 4] )
               frct[m, i, j, 3] = frct[m, i, j, 3]-
                  dsspm * (- 4.0e+00 * rsd[m, i, j, 2]+
                             6.0e+00 * rsd[m, i, j, 3]-
                             4.0e+00 * rsd[m, i, j, 4]+
                                       rsd[m, i, j, 5] )
            end
         end
      end

      for k = 4:nz - 3
         for j = jst:jend
            for i = ist:iend
               for m = 1:5
                  frct[m, i, j, k] = frct[m, i, j, k]-
                     dsspm * (           rsd[m, i, j, k-2]-
                               4.0e+00 * rsd[m, i, j, k-1]+
                               6.0e+00 * rsd[m, i, j, k]-
                               4.0e+00 * rsd[m, i, j, k+1]+
                                         rsd[m, i, j, k+2] )
               end
            end
         end
      end

      for j = jst:jend
         for i = ist:iend
            for m = 1:5
               frct[m, i, j, nz-2] = frct[m, i, j, nz-2]-
                  dsspm * (            rsd[m, i, j, nz-4]-
                             4.0e+00 * rsd[m, i, j, nz-3]+
                             6.0e+00 * rsd[m, i, j, nz-2]-
                             4.0e+00 * rsd[m, i, j, nz-1]  )
               frct[m, i, j, nz-1] = frct[m, i, j, nz-1]-
                  dsspm * (             rsd[m, i, j, nz-3]-
                              4.0e+00 * rsd[m, i, j, nz-2]+
                              5.0e+00 * rsd[m, i, j, nz-1]  )
            end
         end
      end

      return nothing
end
