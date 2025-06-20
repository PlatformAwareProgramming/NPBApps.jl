#---------------------------------------------------------------------
#   compute the right hand side based on exact solution
#---------------------------------------------------------------------

 function exchange_6(g, jbeg, jfin1, isiz03, ny, nz, west, east, comm_solve)

      dum = Array{Float64}(undef, isiz03)

      msgid3 = Ref{MPI.Request}()

#---------------------------------------------------------------------
#   communicate in the east and west directions
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#   receive from east
#---------------------------------------------------------------------
      if jfin1 == ny

        msgid3[] = MPI.Irecv!(view(dum, 1:nz), east, from_e, comm_solve)

        MPI.Wait(msgid3[])

        @inbounds for k = 1:nz
          g[ny+1,k] = dum[k]
        end

      end

#---------------------------------------------------------------------
#   send west
#---------------------------------------------------------------------
      if jbeg == 1
        @inbounds for k = 1:nz
          dum[k] = g[1,k]
        end

        MPI.Send(view(dum, 1:nz), west, from_e, comm_solve)

      end

      return nothing
end
