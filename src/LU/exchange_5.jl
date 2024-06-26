#---------------------------------------------------------------------
#   compute the right hand side based on exact solution
#---------------------------------------------------------------------

 function exchange_5(g, ibeg, ifin1, isiz03, nx, nz, north, south, comm_solve)

      msgid1 = Ref{MPI.Request}()

      dum = Array{Float64}(undef, isiz03)

#---------------------------------------------------------------------
#   communicate in the south and north directions
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#   receive from south
#---------------------------------------------------------------------
      if ifin1 == nx
        msgid1[] = MPI.Irecv!(view(dum, 1:nz), south, from_s, comm_solve)

        MPI.Wait(msgid1[])

        for k = 1:nz
          g[nx+1,k] = dum[k]
        end

      end

#---------------------------------------------------------------------
#   send north
#---------------------------------------------------------------------
      if ibeg == 1
        for k = 1:nz
          dum[k] = g[1,k]
        end

        MPI.Send(view(dum, 1:nz), north, from_s, comm_solve)

      end

      return nothing
end
