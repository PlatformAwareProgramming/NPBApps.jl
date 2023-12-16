using FortranFiles
using OffsetArrays
using Parameters
using Printf


#---------------------------------------------------------------------
#---------------------------------------------------------------------

      function exchange_4(g, h, ibeg, ifin1, jbeg, jfin1)

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
#      integer ibeg, ifin1, jbeg, jfin1
#      DOUBLEPRECISION  g[0:isiz2+1,0:isiz3+1],  
#              h[0:isiz2+1,0:isiz3+1]

#---------------------------------------------------------------------
#  local variables
#---------------------------------------------------------------------
#      integer i, j
#      integer ny2
#      DOUBLEPRECISION  dum[2*isiz02+4]

#      integer msgid1, msgid3
#      integer STATUS[MPI_STATUS_SIZE]
#      integer IERROR



      ny2 = ny + 2

#---------------------------------------------------------------------
#   communicate in the east and west directions
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#   receive from east
#---------------------------------------------------------------------
      if jfin1 == ny
        MPI.Irecv!( dum,
                        2*nx,
                        dp_type,
                        east,
                        from_e,
                        comm_solve,
                        msgid3,
                        IERROR )

        MPI_WAIT( msgid3, STATUS, IERROR )

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

        MPI.Send( dum,
                       2*nx,
                       dp_type,
                       west,
                       from_e,
                       comm_solve,
                       IERROR )

      end

#---------------------------------------------------------------------
#   communicate in the south and north directions
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#   receive from south
#---------------------------------------------------------------------
      if ifin1 == nx
        MPI.Irecv!( dum,
                        2*ny2,
                        dp_type,
                        south,
                        from_s,
                        comm_solve,
                        msgid1,
                        IERROR )

        MPI_WAIT( msgid1, STATUS, IERROR )

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

        MPI.Send( dum,
                       2*ny2,
                       dp_type,
                       north,
                       from_s,
                       comm_solve,
                       IERROR )

      end

      return nothing
      end
