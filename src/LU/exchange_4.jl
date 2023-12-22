#---------------------------------------------------------------------
#   compute the right hand side based on exact solution
#---------------------------------------------------------------------

 function exchange_4(g, h, ibeg, ifin1, jbeg, jfin1)

      dum = Array{Float64}(undef, 2*isiz02+4)

      msgid1 = Ref{MPI.Request}()
      msgid3 = Ref{MPI.Request}()

      ny2 = ny + 2

#---------------------------------------------------------------------
#   communicate in the east and west directions
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#   receive from east
#---------------------------------------------------------------------
      if jfin1 == ny
        msgid3[] = MPI.Irecv!(view(dum, 1:2*nx), east, from_e, comm_solve)

        MPI.Wait(msgid3[])

        for i = 1:nx
          g[i,ny+1] = dum[i]
          h[i,ny+1] = dum[i+nx]
        end

      end

#---------------------------------------------------------------------
#   send west
#---------------------------------------------------------------------
      if jbeg == 1
        for i = 1:nx
          dum[i] = g[i,1]
          dum[i+nx] = h[i,1]
        end

        MPI.Send(view(dum, 1:2*nx), west, from_e, comm_solve)

      end

#---------------------------------------------------------------------
#   communicate in the south and north directions
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#   receive from south
#---------------------------------------------------------------------
      if ifin1 == nx

        msgid1[] = MPI.Irecv!(view(dum, 1:2*ny2), south, from_s, comm_solve)

        MPI.Wait(msgid1[])

        for j = 0:ny+1
          g[nx+1,j] = dum[j+1]
          h[nx+1,j] = dum[j+ny2+1]
        end

      end

#---------------------------------------------------------------------
#   send north
#---------------------------------------------------------------------
      if ibeg == 1
        for j = 0:ny+1
          dum[j+1] = g[1,j]
          dum[j+ny2+1] = h[1,j]
        end

        MPI.Send(view(dum, 1:2*ny2), north, from_s, comm_solve)

      end

      return nothing
end
