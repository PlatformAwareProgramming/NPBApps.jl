#---------------------------------------------------------------------
#---------------------------------------------------------------------

 function exchange_1( g, k, iex,
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

      if iex == 0

          if north != -1
              MPI.Recv!(view(buf1, 1:5, jst:jend), north, from_n, comm_solve)
              @inbounds for j = jst:jend
                  g[1,0,j,k] = buf1[1, j]
                  g[2,0,j,k] = buf1[2, j]
                  g[3,0,j,k] = buf1[3, j]
                  g[4,0,j,k] = buf1[4, j]
                  g[5,0,j,k] = buf1[5, j]
              end
          end

          if west != -1
              MPI.Recv!(view(buf1, 1:5, ist:iend), west, from_w, comm_solve)
              @inbounds for i = ist:iend
                  g[1,i,0,k] = buf1[1, i]
                  g[2,i,0,k] = buf1[2, i]
                  g[3,i,0,k] = buf1[3, i]
                  g[4,i,0,k] = buf1[4, i]
                  g[5,i,0,k] = buf1[5, i]
              end
          end

      elseif iex == 1

          if south != -1
              MPI.Recv!(view(buf1, 1:5, jst:jend), south, from_s, comm_solve)
              @inbounds for j = jst:jend
                  g[1,nx+1,j,k] = buf1[1, j]
                  g[2,nx+1,j,k] = buf1[2, j]
                  g[3,nx+1,j,k] = buf1[3, j]
                  g[4,nx+1,j,k] = buf1[4, j]
                  g[5,nx+1,j,k] = buf1[5, j]
              end
          end

          if east != -1
              MPI.Recv!(view(buf1, 1:5, ist:iend), east, from_e,  comm_solve)
              @inbounds for i = ist:iend
                  g[1,i,ny+1,k] = buf1[1, i]
                  g[2,i,ny+1,k] = buf1[2, i]
                  g[3,i,ny+1,k] = buf1[3, i]
                  g[4,i,ny+1,k] = buf1[4, i]
                  g[5,i,ny+1,k] = buf1[5, i]
              end
          end

      elseif iex == 2

          if south != -1
              @inbounds for j = jst:jend
                  buf[1, j] = g[1,nx,j,k]
                  buf[2, j] = g[2,nx,j,k]
                  buf[3, j] = g[3,nx,j,k]
                  buf[4, j] = g[4,nx,j,k]
                  buf[5, j] = g[5,nx,j,k]
             end
             MPI.Send(view(buf, 1:5, jst:jend), south, from_n, comm_solve)
          end

          if east != -1
             @inbounds for i = ist:iend
                  buf[1, i] = g[1,i,ny,k]
                  buf[2, i] = g[2,i,ny,k]
                  buf[3, i] = g[3,i,ny,k]
                  buf[4, i] = g[4,i,ny,k]
                  buf[5, i] = g[5,i,ny,k]
              end
              MPI.Send(view(buf, 1:5, ist:iend), east, from_w, comm_solve)
          end

      else

          if north != -1
              @inbounds for j = jst:jend
                  buf[1, j] = g[1,1,j,k]
                  buf[2, j] = g[2,1,j,k]
                  buf[3, j] = g[3,1,j,k]
                  buf[4, j] = g[4,1,j,k]
                  buf[5, j] = g[5,1,j,k]
              end
              MPI.Send(view(buf, 1:5, jst:jend), north, from_s, comm_solve)
          end

          if west != -1
              @inbounds for i = ist:iend
                  buf[1, i] = g[1,i,1,k]
                  buf[2, i] = g[2,i,1,k]
                  buf[3, i] = g[3,i,1,k]
                  buf[4, i] = g[4,i,1,k]
                  buf[5, i] = g[5,i,1,k]
              end
              MPI.Send(view(buf, 1:5, ist:iend), west, from_e, comm_solve)
          end

      end

      return nothing
end