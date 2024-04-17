#---------------------------------------------------------------------
#   compute the right hand side based on exact solution
#---------------------------------------------------------------------

 function exchange_5(g, ibeg, ifin1, nz0, nx, nz, south, north, comm_solve)

      msgid1 = Ref{MPI.Request}(MPI.REQUEST_NULL)

      dum = Array{Float64}(undef, nz0)

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
