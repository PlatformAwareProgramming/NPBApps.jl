#---------------------------------------------------------------------
#   to perform pseudo-time stepping SSOR iterations
#   for five nonlinear pde's.
#---------------------------------------------------------------------

const delunm = Array{Float64}(undef, 5) 
const tmat = Array{Float64}(undef, 5, 5)

 function ssor(comm_solve, 
               u,
               rsd,
               frct,
               flux,
               a,
               b,
               c,
               d,
               buf,
               buf1,
               south,
               east,
               north,
               west,
               dt,
               omega,
               nx0,
               ny0,
               nz0,
               timeron,
               tx1,
               tx2,
               tx3,
               ty1,
               ty2,
               ty3,
               tz1,
               tz2,
               tz3,
               nx,
               ny,
               nz,
               ist,
               iend,
               jst,
               jend,
               )

      tv = Array{Float64}(undef, 5, nx0) 

      #@info "$clusterid/$node: SSOR 1"

#---------------------------------------------------------------------
#   begin pseudo-time stepping iterations
#---------------------------------------------------------------------
      tmp = 1.0e+00 / ( omega * ( 2.0e+00 - omega ) )

#---------------------------------------------------------------------
#   perform SSOR iteration
#---------------------------------------------------------------------
         for k = 2:nz - 1
            for j = jst:jend
               for i = ist:iend
                  for m = 1:5
                     rsd[m, i, j, k] = dt * rsd[m, i, j, k]
                  end
               end
            end
         end

         #@info "$clusterid/$node: SSOR 2"

         for k = 2:nz -1
#---------------------------------------------------------------------
#   receive data from north and west
#---------------------------------------------------------------------
            if (timeron) timer_start(t_lcomm) end
            iex = 0
            exchange_1( rsd, k, iex,
                        comm_solve, 
                        buf,
                        buf1,
                        south,
                        east,
                        north,
                        west,
                        nx,
                        ny,
                        ist,
                        iend,
                        jst,
                        jend, 
                     )
            if (timeron) timer_stop(t_lcomm) end

            #@info "$clusterid/$node: SSOR 3 k=$k"

            if (timeron) timer_start(t_blts) end
            for j = jst:jend

#---------------------------------------------------------------------
#   form the lower triangular part of the jacobian matrix
#---------------------------------------------------------------------
             jacld(j, k,
                     u,
                     a,
                     b,
                     c,
                     d,
                     dt,
                     tx1,
                     tx2,
                     ty1,
                     ty2,
                     tz1,
                     tz2,
                     ist,
                     iend,
                    )

            #@info "$clusterid/$node: SSOR 4 k=$k"

#---------------------------------------------------------------------
#   perform the lower triangular solution
#---------------------------------------------------------------------
             blts( j, k, omega, rsd, a, b, c, d, ist, iend)

             #@info "$clusterid/$node: SSOR 5.1 k=$k"
      
            end
            if (timeron) timer_stop(t_blts) end

#---------------------------------------------------------------------
#   send data to east and south
#---------------------------------------------------------------------
            if (timeron) timer_start(t_lcomm) end
            iex = 2
            exchange_1( rsd, k, iex,
                        comm_solve, 
                        buf,
                        buf1,
                        south,
                        east,
                        north,
                        west,
                        nx,
                        ny,
                        ist,
                        iend,
                        jst,
                        jend,
                     )
            if (timeron) timer_stop(t_lcomm) end
            #@info "$clusterid/$node: SSOR 5.2 k=$k"
         end


         for k = nz - 1:-1:2
            #@info "$clusterid/$node: SSOR 6 k=$k"
#---------------------------------------------------------------------
#   receive data from south and east
#---------------------------------------------------------------------
            if (timeron) timer_start(t_ucomm) end
            iex = 1
            exchange_1( rsd, k, iex,
                        comm_solve, 
                        buf,
                        buf1,
                        south,
                        east,
                        north,
                        west,
                        nx,
                        ny,
                        ist,
                        iend,
                        jst,
                        jend,
                       )
            if (timeron) timer_stop(t_ucomm) end

            #@info "$clusterid/$node: SSOR 7 k=$k"

            if (timeron) timer_start(t_buts) end
            for j = jend:-1:jst

#---------------------------------------------------------------------
#   form the strictly upper triangular part of the jacobian matrix
#---------------------------------------------------------------------
              jacu(j, k, 
                     u,
                     a,
                     b,
                     c,
                     d,
                     dt,
                     tx1,
                     tx2,
                     ty1,
                     ty2,
                     tz1,
                     tz2,
                     ist,
                     iend,
                   )

               #@info "$clusterid/$node: SSOR 8 k=$k"
#---------------------------------------------------------------------
#   perform the upper triangular solution
#---------------------------------------------------------------------
               buts( j, k,
                          omega,
                          rsd, tv,
                          d, a, b, c,
                          ist, iend)

               #@info "$clusterid/$node: SSOR 9 k=$k"                          
            end
            if (timeron) timer_stop(t_buts) end

#---------------------------------------------------------------------
#   send data to north and west
#---------------------------------------------------------------------
            if (timeron) timer_start(t_ucomm) end
            iex = 3
            exchange_1( rsd, k, iex,
                        comm_solve, 
                        buf,
                        buf1,
                        south,
                        east,
                        north,
                        west,
                        nx,
                        ny,
                        ist,
                        iend,
                        jst,
                        jend, 
                     )
            if (timeron) timer_stop(t_ucomm) end
            #@info "$clusterid/$node: SSOR 10 k=$k"
         end

#---------------------------------------------------------------------
#   update the variables
#---------------------------------------------------------------------

         for k = 2:nz-1
            for j = jst:jend
               for i = ist:iend
                  for m = 1:5
                     u[ m, i, j, k ] = u[ m, i, j, k ] + tmp * rsd[ m, i, j, k ]
                  end
               end
            end
         end

         #@info "$clusterid/$node: SSOR 12"

#---------------------------------------------------------------------
#   compute the steady-state residuals
#---------------------------------------------------------------------
         rhs(               
               comm_solve, 
               u,
               rsd,
               frct,
               flux,
               buf,
               buf1,
               south,
               east,
               north,
               west,
               timeron,
               tx1,
               tx2,
               tx3,
               ty1,
               ty2,
               ty3,
               tz1,
               tz2,
               tz3,
               nx,
               ny,
               nz,
               ist,
               iend,
               jst,
               jend,
            )

      return nothing
end
