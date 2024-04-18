#---------------------------------------------------------------------
#---------------------------------------------------------------------

 function pintgr(u, phi1, phi2, ny0, nz0, nx, ny, nz, ipt, jpt, west, east, south, north, comm_solve)

     @info "$clusterid/$node: pintgr 1"
#---------------------------------------------------------------------
#   set up the sub-domains for integration in each processor
#---------------------------------------------------------------------
      ibeg = nx + 1
      ifin = 0
      iglob1 = ipt + 1
      iglob2 = ipt + nx
      if (iglob1 >= ii1 && iglob2 < ii2+nx) ibeg = 1 end
      if (iglob1 > ii1-nx && iglob2 <= ii2) ifin = nx end
      if (ii1 >= iglob1 && ii1 <= iglob2) ibeg = ii1 - ipt end
      if (ii2 >= iglob1 && ii2 <= iglob2) ifin = ii2 - ipt end
      jbeg = ny + 1
      jfin = 0
      jglob1 = jpt + 1
      jglob2 = jpt + ny
      if (jglob1 >= ji1 && jglob2 < ji2+ny) jbeg = 1 end
      if (jglob1 > ji1-ny && jglob2 <= ji2) jfin = ny end
      if (ji1 >= jglob1 && ji1 <= jglob2) jbeg = ji1 - jpt end
      if (ji2 >= jglob1 && ji2 <= jglob2) jfin = ji2 - jpt end
      ifin1 = ifin
      jfin1 = jfin
      if (ipt + ifin1 == ii2) ifin1 = ifin -1 end
      if (jpt + jfin1 == ji2) jfin1 = jfin -1 end

      @info "$clusterid/$node: pintgr 2"
      
#---------------------------------------------------------------------
#   initialize
#---------------------------------------------------------------------
      for k = 0:nz+1
        for i = 0:ny+1
          phi1[i, k] = 0.
          phi2[i, k] = 0.
        end
      end

      @info "$clusterid/$node: pintgr 3   ny0=$ny0, nx=$nx, ny=$ny, ibeg=$ibeg, ifin1=$ifin1, jbeg=$jbeg, jfin1=$jfin1"
      
      for j = jbeg:jfin
         jglob = jpt + j
         for i = ibeg:ifin
            iglob = ipt + i

            k = ki1
            phi1[i, j] = c2*(  u[5, i, j, k]-
                  0.50e+00 * (  u[2, i, j, k] ^ 2+
                                u[3, i, j, k] ^ 2+
                                u[4, i, j, k] ^ 2 )/
                               u[1, i, j, k] )

            k = ki2

            phi2[i, j] = c2*(  u[5, i, j, k]-
                  0.50e+00 * (  u[2, i, j, k] ^ 2+
                                u[3, i, j, k] ^ 2+
                                u[4, i, j, k] ^ 2 )/
                               u[1, i, j, k] )
         end
      end

      @info "$clusterid/$node: pintgr 4"

#---------------------------------------------------------------------
#  communicate in i and j directions
#---------------------------------------------------------------------
      exchange_4(phi1, phi2, ibeg, ifin1, jbeg, jfin1, ny0, nx, ny, west, east, south, north, comm_solve)

      @info "$clusterid/$node: pintgr 5"
      
      frc1 = 0.0e+00

      for j = jbeg:jfin1
         for i = ibeg:ifin1
            frc1 += (  phi1[i, j]+
                            phi1[i+1, j]+
                            phi1[i, j+1]+
                            phi1[i+1, j+1]+
                            phi2[i, j]+
                            phi2[i+1, j]+
                            phi2[i, j+1]+
                            phi2[i+1, j+1] )
         end
      end

      @info "$clusterid/$node: pintgr 6"
      
#---------------------------------------------------------------------
#  compute the global sum of individual contributions to frc1
#---------------------------------------------------------------------
      #dummy = frc1
      #MPI_ALLREDUCE( dummy, frc1, 1, dp_type, MPI_SUM, comm_solve, IERROR )

      frc1 = MPI.Allreduce(frc1, MPI.SUM, comm_solve)

      @info "$clusterid/$node: pintgr 7"
      
      frc1 = dxi * deta * frc1

#---------------------------------------------------------------------
#   initialize
#---------------------------------------------------------------------
      for k = 0:nz+1
        for i = 0:ny+1
          phi1[i, k] = 0.
          phi2[i, k] = 0.
        end
      end
      jglob = jpt + jbeg
      ind1 = 0
      if jglob == ji1
        ind1 = 1
        for k = ki1:ki2
           for i = ibeg:ifin
              iglob = ipt + i
              phi1[i, k] = c2*(  u[5, i, jbeg, k]-
                    0.50e+00 * (  u[2, i, jbeg, k] ^ 2+
                                  u[3, i, jbeg, k] ^ 2+
                                  u[4, i, jbeg, k] ^ 2 )/
                                 u[1, i, jbeg, k] )
           end
        end
      end

      @info "$clusterid/$node: pintgr 8"
      
      jglob = jpt + jfin
      ind2 = 0
      if jglob == ji2
        ind2 = 1
        for k = ki1:ki2
           for i = ibeg:ifin
              iglob = ipt + i
              phi2[i, k] = c2*(  u[5, i, jfin, k]-
                    0.50e+00 * (  u[2, i, jfin, k] ^ 2+
                                  u[3, i, jfin, k] ^ 2+
                                  u[4, i, jfin, k] ^ 2 )/
                                 u[1, i, jfin, k] )
           end
        end
      end

      @info "$clusterid/$node: pintgr 9"
      
#---------------------------------------------------------------------
#  communicate in i direction
#---------------------------------------------------------------------
      if ind1 == 1
        exchange_5(phi1, ibeg, ifin1, nz0, nx, nz, south, north, comm_solve)
      end

      @info "$clusterid/$node: pintgr 10"
      
      if ind2 == 1
        exchange_5(phi2, ibeg, ifin1, nz0, nx, nz, south, north, comm_solve)
      end

      @info "$clusterid/$node: pintgr 11"
      
      frc2 = 0.0e+00
      for k = ki1:ki2-1
         for i = ibeg:ifin1
            frc2 += (  phi1[i, k]+
                            phi1[i+1, k]+
                            phi1[i, k+1]+
                            phi1[i+1, k+1]+
                            phi2[i, k]+
                            phi2[i+1, k]+
                            phi2[i, k+1]+
                            phi2[i+1, k+1] )
         end
      end

      @info "$clusterid/$node: pintgr 12"
      
#---------------------------------------------------------------------
#  compute the global sum of individual contributions to frc2
#---------------------------------------------------------------------
      #dummy = frc2
      #MPI_ALLREDUCE( dummy, frc2, 1, dp_type, MPI_SUM, comm_solve, IERROR )
      
      frc2 = MPI.Allreduce(frc2, MPI.SUM, comm_solve)

      @info "$clusterid/$node: pintgr 13"
      
      frc2 = dxi * dzeta * frc2

#---------------------------------------------------------------------
#   initialize
#---------------------------------------------------------------------
      for k = 0:nz+1
        for i = 0:ny+1
          phi1[i, k] = 0.
          phi2[i, k] = 0.
        end
      end
      iglob = ipt + ibeg
      ind1 = 0
      if iglob == ii1
        ind1 = 1
        for k = ki1:ki2
           for j = jbeg:jfin
              jglob = jpt + j
              phi1[j, k] = c2*(  u[5, ibeg, j, k]-
                    0.50e+00 * (  u[2, ibeg, j, k] ^ 2+
                                  u[3, ibeg, j, k] ^ 2+
                                  u[4, ibeg, j, k] ^ 2 )/
                                 u[1, ibeg, j, k] )
           end
        end
      end

      @info "$clusterid/$node: pintgr 14"
      
      iglob = ipt + ifin
      ind2 = 0
      if iglob == ii2
        ind2 = 1
        for k = ki1:ki2
           for j = jbeg:jfin
              jglob = jpt + j
              phi2[j, k] = c2*(  u[5, ifin, j, k]-
                    0.50e+00 * (  u[2, ifin, j, k] ^ 2+
                                  u[3, ifin, j, k] ^ 2+
                                  u[4, ifin, j, k] ^ 2 )/
                                 u[1, ifin, j, k] )
           end
        end
      end

      @info "$clusterid/$node: pintgr 15"
      
#---------------------------------------------------------------------
#  communicate in j direction
#---------------------------------------------------------------------
      if ind1 == 1
        exchange_6(phi1, jbeg, jfin1, nz0, ny, nz, east, west, comm_solve)
      end
     
      @info "$clusterid/$node: pintgr 16"
      
      if ind2 == 1
        exchange_6(phi2, jbeg, jfin1, nz0, ny, nz, east, west, comm_solve)
      end

      @info "$clusterid/$node: pintgr 17"
      
      frc3 = 0.0e+00

      for k = ki1:ki2-1
         for j = jbeg:jfin1
            frc3 += (  phi1[j, k]+
                            phi1[j+1, k]+
                            phi1[j, k+1]+
                            phi1[j+1, k+1]+
                            phi2[j, k]+
                            phi2[j+1, k]+
                            phi2[j, k+1]+
                            phi2[j+1, k+1] )
         end
      end

      @info "$clusterid/$node: pintgr 18"
      
#---------------------------------------------------------------------
#  compute the global sum of individual contributions to frc3
#---------------------------------------------------------------------
      #dummy = frc3
      #MPI_ALLREDUCE( dummy, frc3, 1, dp_type, MPI_SUM, comm_solve, IERROR )

      frc3 = MPI.Allreduce(frc3, MPI.SUM, comm_solve)

      @info "$clusterid/$node: pintgr 19"
      
      frc3 = deta * dzeta * frc3
      frc = 0.25e+00 * ( frc1 + frc2 + frc3 )
#      if (id.eq.0) write (*,1001) frc

      @info "$clusterid/$node: pintgr 20"
            
      return frc
end
