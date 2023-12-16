using FortranFiles
using OffsetArrays
using Parameters
using Printf


#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function exchange_6(g, jbeg, jfin1)

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
#      integer jbeg, jfin1
#      DOUBLEPRECISION  g[0:isiz2+1,0:isiz3+1]

#---------------------------------------------------------------------
#  local parameters
#---------------------------------------------------------------------
#      integer k
#      DOUBLEPRECISION  dum[isiz03]

#      integer msgid3
#      integer STATUS[MPI_STATUS_SIZE]
#      integer IERROR



#---------------------------------------------------------------------
#   communicate in the east and west directions
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#   receive from east
#---------------------------------------------------------------------
      if jfin1 == ny
        MPI.Irecv!( dum,
                        nz,
                        dp_type,
                        east,
                        from_e,
                        comm_solve,
                        msgid3,
                        IERROR )

        MPI_WAIT( msgid3, STATUS, IERROR )

        for k = 1:nz
          g[ny+1,k] = dum[k]
        end

      end

#---------------------------------------------------------------------
#   send west
#---------------------------------------------------------------------
      if jbeg == 1
        for k = 1:nz
          dum[k] = g[1,k]
        end

        MPI.Send( dum,
                       nz,
                       dp_type,
                       west,
                       from_e,
                       comm_solve,
                       IERROR )

      end

      return nothing
      end
