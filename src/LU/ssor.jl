#---------------------------------------------------------------------
#   to perform pseudo-time stepping SSOR iterations
#   for five nonlinear pde's.
#---------------------------------------------------------------------

const delunm = Array{Float64}(undef, 5) 
const tmat = Array{Float64}(undef, 5, 5)

 function ssor(niter,
               comm_solve, 
               id,
               rsdnm,
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
               isiz1,
               isiz2,
               isiz3,
               inorm,
               itmax,
               dt,
               omega,
               tolrsd,
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
               ipt,
               ny,
               jpt,
               nz,
               ist,
               iend,
               jst,
               jend,
               )

      tv = Array{Float64}(undef, 5, isiz1) 

#---------------------------------------------------------------------
#   begin pseudo-time stepping iterations
#---------------------------------------------------------------------
      tmp = 1.0e+00 / ( omega * ( 2.0e+00 - omega ) )

#---------------------------------------------------------------------
#   initialize a,b,c,d to zero (guarantees that page tables have been
#   formed, if applicable on given architecture, before timestepping).
#---------------------------------------------------------------------
      for i = 1:isiz1
         for m = 1:5
            for k = 1:5
               a[k, m, i] = 0.0e0
               b[k, m, i] = 0.0e0
               c[k, m, i] = 0.0e0
               d[k, m, i] = 0.0e0
            end
         end
      end

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

#---------------------------------------------------------------------
#   compute the L2 norms of newton iteration residuals
#---------------------------------------------------------------------
      l2norm( isiz1, isiz2, isiz3, nx0, ny0, nz0,
                   ist, iend, jst, jend,
                   rsd, rsdnm, 
                   comm_solve, 
                   timeron,
            )

      for i = 1:t_last
         timer_clear(i)
      end

      MPI.Barrier(comm_solve)

      timer_clear(1)
      timer_start(1)

#---------------------------------------------------------------------
#   the timestep loop
#---------------------------------------------------------------------
      for istep = 1:niter

         if id == 0
            if mod(istep, 20) == 0 || istep == itmax || istep == 1
               if (niter > 1) 
                  @printf(stdout, " Time step %4i\n", istep) 
               end
            end
         end

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

#---------------------------------------------------------------------
#   perform the lower triangular solution
#---------------------------------------------------------------------
             blts( isiz1, isiz2, isiz3,
                          nx, ny, nz, j, k,
                          omega,
                          rsd,
                          a, b, c, d,
                          ist, iend, jst, jend,
                          nx0, ny0, ipt, jpt, tmat)
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
         end

         for k = nz - 1:-1:2
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

#---------------------------------------------------------------------
#   perform the upper triangular solution
#---------------------------------------------------------------------
               buts( isiz1, isiz2, isiz3,
                          nx, ny, nz, j, k,
                          omega,
                          rsd, tv,
                          d, a, b, c,
                          ist, iend, jst, jend,
                          nx0, ny0, ipt, jpt, tmat)
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

#---------------------------------------------------------------------
#   compute the max-norms of newton iteration corrections
#---------------------------------------------------------------------
         if mod(istep, inorm) == 0
            l2norm(isiz1, isiz2, isiz3, nx0, ny0, nz0, ist, iend, jst, jend, rsd, delunm, 
                   comm_solve, 
                   timeron,
                  )
#            if ( ipr .eq. 1 .and. id .eq. 0 ) then
#                write (*,1006) ( delunm(m), m = 1, 5 )
#            else if ( ipr .eq. 2 .and. id .eq. 0 ) then
#                write (*,'(i5,f15.6)') istep,delunm(5)
#            end if
         end

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

#---------------------------------------------------------------------
#   compute the max-norms of newton iteration residuals
#---------------------------------------------------------------------
         if (mod(istep, inorm) == 0) || (istep == itmax)
            l2norm( isiz1, isiz2, isiz3, nx0, ny0, nz0, ist, iend, jst, jend, rsd, rsdnm, 
                    comm_solve, 
                    timeron,
                  )
#            if ( ipr .eq. 1.and.id.eq.0 ) then
#                write (*,1007) ( rsdnm[m], m = 1, 5 )
#            end if
         end

#---------------------------------------------------------------------
#   check the newton-iteration residuals against the tolerance levels
#---------------------------------------------------------------------
         if ( rsdnm[1] < tolrsd[1] ) &&(
               rsdnm[2] < tolrsd[2] ) &&(
               rsdnm[3] < tolrsd[3] ) &&(
               rsdnm[4] < tolrsd[4] ) &&(
               rsdnm[5] < tolrsd[5] )
            if id == 0
               @printf(stdout, " \n convergence was achieved after %4i pseudo-time steps\n", istep)
            end
             @goto L900
         end

      end
      @label L900

      timer_stop(1)
      wtime = timer_read(1)


      #MPI_ALLREDUCE( wtime, maxtime, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm_solve, IERROR )
      
      global maxtime = MPI.Allreduce(wtime, MPI.MAX, comm_solve)

# 1001 format(1x/5x,'pseudo-time SSOR iteration no.=',i4/)
# 1004 format(1x/1x,'convergence was achieved after ',i4,  		         ' pseudo-time steps' )
# 1006 format(1x/1x,'RMS-norm of SSOR-iteration correction ',  		       'for first pde  = ',1pe12.5/,  		       1x,'RMS-norm of SSOR-iteration correction ',  		       'for second pde = ',1pe12.5/,  		       1x,'RMS-norm of SSOR-iteration correction ',  		       'for third pde  = ',1pe12.5/,  		       1x,'RMS-norm of SSOR-iteration correction ',  		       'for fourth pde = ',1pe12.5/,  		       1x,'RMS-norm of SSOR-iteration correction ',  		       'for fifth pde  = ',1pe12.5)
# 1007 format(1x/1x,'RMS-norm of steady-state residual for ',  		       'first pde  = ',1pe12.5/,  		       1x,'RMS-norm of steady-state residual for ',  		       'second pde = ',1pe12.5/,  		       1x,'RMS-norm of steady-state residual for ',  		       'third pde  = ',1pe12.5/,  		       1x,'RMS-norm of steady-state residual for ',  		       'fourth pde = ',1pe12.5/,  		       1x,'RMS-norm of steady-state residual for ',  		       'fifth pde  = ',1pe12.5)

      return nothing
end
