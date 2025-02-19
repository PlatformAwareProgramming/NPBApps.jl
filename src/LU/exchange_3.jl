#---------------------------------------------------------------------
#   compute the right hand side based on exact solution
#---------------------------------------------------------------------

const mid = Ref{MPI.Request}(MPI.REQUEST_NULL)

function exchange_3(g, iex,
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

      if iex == 0
        
#---------------------------------------------------------------------
#   communicate in the south and north directions
#---------------------------------------------------------------------
        if north != -1
            mid[] = MPI.Irecv!(view(buf1, 1:10*ny*nz), north, from_n, comm_solve)
        end

#---------------------------------------------------------------------
#   send south
#---------------------------------------------------------------------
        if south != -1
            for k = 1:nz
              for j = 1:ny
                ipos1 = (k-1)*ny + j
                ipos2 = ipos1 + ny*nz
                buf[1, ipos1] = g[1,nx-1,j,k]
                buf[2, ipos1] = g[2,nx-1,j,k]
                buf[3, ipos1] = g[3,nx-1,j,k]
                buf[4, ipos1] = g[4,nx-1,j,k]
                buf[5, ipos1] = g[5,nx-1,j,k]
                buf[1, ipos2] = g[1,nx,j,k]
                buf[2, ipos2] = g[2,nx,j,k]
                buf[3, ipos2] = g[3,nx,j,k]
                buf[4, ipos2] = g[4,nx,j,k]
                buf[5, ipos2] = g[5,nx,j,k]
              end
            end

            MPI.Send(view(buf, 1:10*ny*nz), south, from_n, comm_solve)
        end

#---------------------------------------------------------------------
#   receive from north
#---------------------------------------------------------------------
        if north != -1
          MPI.Wait(mid[])

          for k = 1:nz
            for j = 1:ny
              ipos1 = (k-1)*ny + j
              ipos2 = ipos1 + ny*nz
              g[1,-1,j,k] = buf1[1, ipos1]
              g[2,-1,j,k] = buf1[2, ipos1]
              g[3,-1,j,k] = buf1[3, ipos1]
              g[4,-1,j,k] = buf1[4, ipos1]
              g[5,-1,j,k] = buf1[5, ipos1]
              g[1,0,j,k] = buf1[1, ipos2]
              g[2,0,j,k] = buf1[2, ipos2]
              g[3,0,j,k] = buf1[3, ipos2]
              g[4,0,j,k] = buf1[4, ipos2]
              g[5,0,j,k] = buf1[5, ipos2]
            end
          end

        end

        if south != -1
           mid[] = MPI.Irecv!(view(buf1, 1:10*ny*nz), south, from_s, comm_solve)
        end

#---------------------------------------------------------------------
#   send north
#---------------------------------------------------------------------
        if north != -1
          for k = 1:nz
            for j = 1:ny
              ipos1 = (k-1)*ny + j
              ipos2 = ipos1 + ny*nz
              buf[1, ipos1] = g[1,2,j,k]
              buf[2, ipos1] = g[2,2,j,k]
              buf[3, ipos1] = g[3,2,j,k]
              buf[4, ipos1] = g[4,2,j,k]
              buf[5, ipos1] = g[5,2,j,k]
              buf[1, ipos2] = g[1,1,j,k]
              buf[2, ipos2] = g[2,1,j,k]
              buf[3, ipos2] = g[3,1,j,k]
              buf[4, ipos2] = g[4,1,j,k]
              buf[5, ipos2] = g[5,1,j,k]
            end
          end

          MPI.Send(view(buf, 1:10*ny*nz), north, from_s, comm_solve)
        end

#---------------------------------------------------------------------
#   receive from south
#---------------------------------------------------------------------
        if south != -1
          MPI.Wait(mid[])

          for k = 1:nz
            for j = 1:ny
              ipos1 = (k-1)*ny + j
              ipos2 = ipos1 + ny*nz
              g[1,nx+2,j,k] = buf1[1, ipos1]
              g[2,nx+2,j,k] = buf1[2, ipos1]
              g[3,nx+2,j,k] = buf1[3, ipos1]
              g[4,nx+2,j,k] = buf1[4, ipos1]
              g[5,nx+2,j,k] = buf1[5, ipos1]
              g[1,nx+1,j,k] = buf1[1, ipos2]
              g[2,nx+1,j,k] = buf1[2, ipos2]
              g[3,nx+1,j,k] = buf1[3, ipos2]
              g[4,nx+1,j,k] = buf1[4, ipos2]
              g[5,nx+1,j,k] = buf1[5, ipos2]
            end
          end
        end

      else

#---------------------------------------------------------------------
#   communicate in the east and west directions
#---------------------------------------------------------------------
      if west != -1
          mid[] = MPI.Irecv!(view(buf1, 1:10*nx*nz), west, from_w, comm_solve)
      end

#---------------------------------------------------------------------
#   send east
#---------------------------------------------------------------------
        if east != -1
          for k = 1:nz
            for i = 1:nx
              ipos1 = (k-1)*nx + i
              ipos2 = ipos1 + nx*nz
              buf[1, ipos1] = g[1,i,ny-1,k]
              buf[2, ipos1] = g[2,i,ny-1,k]
              buf[3, ipos1] = g[3,i,ny-1,k]
              buf[4, ipos1] = g[4,i,ny-1,k]
              buf[5, ipos1] = g[5,i,ny-1,k]
              buf[1, ipos2] = g[1,i,ny,k]
              buf[2, ipos2] = g[2,i,ny,k]
              buf[3, ipos2] = g[3,i,ny,k]
              buf[4, ipos2] = g[4,i,ny,k]
              buf[5, ipos2] = g[5,i,ny,k]
            end
          end

          MPI.Send(view(buf, 1:10*nx*nz), east, from_w, comm_solve)
        end

#---------------------------------------------------------------------
#   receive from west
#---------------------------------------------------------------------
        if west != -1
          MPI.Wait(mid[])

          for k = 1:nz
            for i = 1:nx
              ipos1 = (k-1)*nx + i
              ipos2 = ipos1 + nx*nz
              g[1,i,-1,k] = buf1[1, ipos1]
              g[2,i,-1,k] = buf1[2, ipos1]
              g[3,i,-1,k] = buf1[3, ipos1]
              g[4,i,-1,k] = buf1[4, ipos1]
              g[5,i,-1,k] = buf1[5, ipos1]
              g[1,i,0,k] = buf1[1, ipos2]
              g[2,i,0,k] = buf1[2, ipos2]
              g[3,i,0,k] = buf1[3, ipos2]
              g[4,i,0,k] = buf1[4, ipos2]
              g[5,i,0,k] = buf1[5, ipos2]
            end
          end

        end

      if east != -1
          mid[] = MPI.Irecv!(view(buf1, 1:10*nx*nz), east, from_e, comm_solve)
      end

#---------------------------------------------------------------------
#   send west
#---------------------------------------------------------------------
      if west != -1
          for k = 1:nz
            for i = 1:nx
              ipos1 = (k-1)*nx + i
              ipos2 = ipos1 + nx*nz
              buf[1, ipos1] = g[1,i,2,k]
              buf[2, ipos1] = g[2,i,2,k]
              buf[3, ipos1] = g[3,i,2,k]
              buf[4, ipos1] = g[4,i,2,k]
              buf[5, ipos1] = g[5,i,2,k]
              buf[1, ipos2] = g[1,i,1,k]
              buf[2, ipos2] = g[2,i,1,k]
              buf[3, ipos2] = g[3,i,1,k]
              buf[4, ipos2] = g[4,i,1,k]
              buf[5, ipos2] = g[5,i,1,k]
            end
          end

          MPI.Send(view(buf, 1:10*nx*nz), west, from_e, comm_solve)
        end

#---------------------------------------------------------------------
#   receive from east
#---------------------------------------------------------------------
        if east != -1
          MPI.Wait(mid[])

          for k = 1:nz
            for i = 1:nx
              ipos1 = (k-1)*nx + i
              ipos2 = ipos1 + nx*nz
              g[1,i,ny+2,k] = buf1[1, ipos1]
              g[2,i,ny+2,k] = buf1[2, ipos1]
              g[3,i,ny+2,k] = buf1[3, ipos1]
              g[4,i,ny+2,k] = buf1[4, ipos1]
              g[5,i,ny+2,k] = buf1[5, ipos1]
              g[1,i,ny+1,k] = buf1[1, ipos2]
              g[2,i,ny+1,k] = buf1[2, ipos2]
              g[3,i,ny+1,k] = buf1[3, ipos2]
              g[4,i,ny+1,k] = buf1[4, ipos2]
              g[5,i,ny+1,k] = buf1[5, ipos2]
            end
          end

        end

      end

      return nothing
end
