using FortranFiles
using OffsetArrays
using Parameters
using Printf


#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function exchange_5(g, ibeg, ifin1)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#   compute the right hand side based on exact solution
#---------------------------------------------------------------------

#      use lu_data
#      use mpinpb

#      implicit none

#---------------------------------------------------------------------
#  input parameters
#---------------------------------------------------------------------
#      integer ibeg, ifin1
#      DOUBLEPRECISION  g[0:isiz2+1,0:isiz3+1]

#---------------------------------------------------------------------
#  local variables
#---------------------------------------------------------------------
#      integer k
#      DOUBLEPRECISION  dum[isiz03]

#      integer msgid1
#      integer STATUS[MPI_STATUS_SIZE]
#      integer IERROR



#---------------------------------------------------------------------
#   communicate in the south and north directions
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#   receive from south
#---------------------------------------------------------------------
      if ifin1 == nx
        MPI.Irecv!( dum,
                        nz,
                        dp_type,
                        south,
                        from_s,
                        comm_solve,
                        msgid1,
                        IERROR )

        MPI_WAIT( msgid1, STATUS, IERROR )

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

        MPI.Send( dum,
                       nz,
                       dp_type,
                       north,
                       from_s,
                       comm_solve,
                       IERROR )

      end

      return nothing
      end
