#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#  bt_data module
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# The following include file is generated automatically by the
# "setparams" utility. It defines 
#      maxcells:      the square root of the maximum number of processors
#      problem_size:  12, 64, 102, 162 (for class S, A, B, C)
#      dt_default:    default time step for this problem size if no
#                     config file
#      niter_default: default number of iterations for this problem size
#---------------------------------------------------------------------

const aa = 1
const bb = 2
const cc = 3
const BLOCK_SIZE = 5

#const grid_points = zeros(Integer, 3)

const EAST = 2000 
const WEST = 3000      
const NORTH = 4000 
const SOUTH = 5000
const BOTTOM = 6000 
const TOP = 7000

const predecessor = zeros(Int64, 3)
const successor = zeros(Int64,3)
const grid_size = zeros(Int64, 3)

#---------------------------------------------------------------------
#     These are used by btio
#---------------------------------------------------------------------
#      integer collbuf_nodes, collbuf_size, iosize, eltext,  
#              combined_btype, fp, idump, record_length, element,  
#              combined_ftype, idump_sub, rd_interval
#      DOUBLEPRECISION SUM[niter_default], xce_sub[5]
#      integer(kind=8) :: iseek

#---------------------------------------------------------------------
#     Timer constants
#---------------------------------------------------------------------

const t_total = 1
const t_io = 2
const t_rhs = 3
const t_xsolve = 4
const t_ysolve = 5
const t_zsolve = 6
const t_bpack = 7
const t_exch = 8
const t_xcomm = 9
const t_ycomm = 10
const t_zcomm = 11
const t_comp = 12
const t_comm = 13
const t_enorm = 12
const t_iov = 13
const t_last = 13

#---------------------------------------------------------------------
# allocate space dynamically for data arrays
#---------------------------------------------------------------------

 function alloc_space(maxcells, problem_size)

       MAX_CELL_DIM = div(problem_size, maxcells)+1

       IMAX = MAX_CELL_DIM
       JMAX = MAX_CELL_DIM
       KMAX = MAX_CELL_DIM

       BUF_SIZE = MAX_CELL_DIM*MAX_CELL_DIM*(maxcells-1)*60+1

       cell_coord = zeros(Int64, 3, maxcells)
       cell_low = zeros(Int64, 3, maxcells)
       cell_high = zeros(Int64, 3, maxcells)
       cell_size = zeros(Int64, 3, maxcells)
       cell_start = zeros(Int64, 3, maxcells)
       cell_end = zeros(Int64, 3, maxcells)
       slice = zeros(Int64, 3, maxcells)

      forcing0 = zeros(Float64, 5, IMAX, JMAX, KMAX, maxcells)
      forcing = OffsetArray(forcing0, 1:5, 0:IMAX-1, 0:JMAX-1, 0:KMAX-1, 1:maxcells)            

      u0 = zeros(Float64, 5, IMAX+4, JMAX+4, KMAX + 4, maxcells)
      u = OffsetArray(u0, 1:5, -2:IMAX+1, -2:JMAX+1, -2:KMAX+1, 1:maxcells)

      rhs0 = zeros(Float64, 5, IMAX+1, JMAX+1, KMAX+1, maxcells)
      rhs = OffsetArray(rhs0, 1:5, -1:IMAX-1, -1:JMAX-1, -1:KMAX-1, 1:maxcells)

      lhsc0 = zeros(Float64, 5, 5, IMAX+1, JMAX+1, KMAX+1, maxcells)
      lhsc = OffsetArray(lhsc0, 1:5, 1:5, -1:IMAX-1, -1:JMAX-1, -1:KMAX-1, 1:maxcells)

      backsub_info0 = zeros(Float64, 5, MAX_CELL_DIM+1, MAX_CELL_DIM+1, maxcells)
      backsub_info = OffsetArray(backsub_info0, 1:5, 0:MAX_CELL_DIM, 0:MAX_CELL_DIM, 1:maxcells)

      in_buffer = Array{Float64}(undef, BUF_SIZE)
      out_buffer = Array{Float64}(undef, BUF_SIZE)

      #cv0 = zeros(Float64, MAX_CELL_DIM + 4)
      #cv = OffsetArray(cv0, -2:MAX_CELL_DIM+1)

      #rhon0 = zeros(Float64, MAX_CELL_DIM + 4)
      #rhon = OffsetArray(rhon0, -2:MAX_CELL_DIM+1)

      #rhos0 = zeros(Float64, MAX_CELL_DIM + 4)
      #rhos = OffsetArray(rhos0, -2:MAX_CELL_DIM+1)
      
      #rhoq0 = zeros(Float64, MAX_CELL_DIM + 4)
      #rhoq = OffsetArray(rhoq0, -2:MAX_CELL_DIM+1)

      cuf0 = zeros(Float64, MAX_CELL_DIM + 4)
      cuf = OffsetArray(cuf0, -2:MAX_CELL_DIM+1)

      q0 = zeros(Float64, MAX_CELL_DIM + 4)
      q = OffsetArray(q0, -2:MAX_CELL_DIM+1)
      
      ue0 = zeros(Float64, MAX_CELL_DIM + 4, 5)
      ue = OffsetArray(ue0, -2:MAX_CELL_DIM+1, 1:5)

      buf0 = zeros(Float64, MAX_CELL_DIM + 4, 5)
      buf = OffsetArray(buf0, -2:MAX_CELL_DIM+1, 1:5)

      fjac0 = zeros(Float64, 5, 5, MAX_CELL_DIM + 4)
      fjac = OffsetArray(fjac0, 1:5, 1:5, -2:MAX_CELL_DIM+1)

      njac0 = zeros(Float64, 5, 5, MAX_CELL_DIM + 4)
      njac = OffsetArray(njac0, 1:5, 1:5, -2:MAX_CELL_DIM+1)

      lhsa0 = zeros(Float64, 5, 5, MAX_CELL_DIM + 2)
      lhsa = OffsetArray(lhsa0, 1:5, 1:5, -1:MAX_CELL_DIM)

      lhsb0 = zeros(Float64, 5, 5, MAX_CELL_DIM + 2)
      lhsb = OffsetArray(lhsb0, 1:5, 1:5, -1:MAX_CELL_DIM)            

      us0 = zeros(Float64, IMAX+2, JMAX+2, KMAX+2, maxcells)
      us = OffsetArray(us0, -1:IMAX, -1:JMAX, -1:KMAX, 1:maxcells)
      
      vs0 = zeros(Float64, IMAX+2, JMAX+2, KMAX+2, maxcells)
      vs = OffsetArray(vs0, -1:IMAX, -1:JMAX, -1:KMAX, 1:maxcells)

      ws0 = zeros(Float64, IMAX+2, JMAX+2, KMAX+2, maxcells)
      ws = OffsetArray(ws0, -1:IMAX, -1:JMAX, -1:KMAX, 1:maxcells)

      qs0 = zeros(Float64, IMAX+2, JMAX+2, KMAX+2, maxcells)
      qs = OffsetArray(qs0, -1:IMAX, -1:JMAX, -1:KMAX, 1:maxcells)
      
      rho_i0 = zeros(Float64, IMAX+2, JMAX+2, KMAX+2, maxcells)
      rho_i = OffsetArray(rho_i0, -1:IMAX, -1:JMAX, -1:KMAX, 1:maxcells)

      square0 = zeros(Float64, IMAX+2, JMAX+2, KMAX+2, maxcells)
      square = OffsetArray(square0, -1:IMAX, -1:JMAX, -1:KMAX, 1:maxcells)
      
      return  MAX_CELL_DIM, IMAX, JMAX, KMAX, BUF_SIZE, 
              cell_coord, cell_low, cell_high, cell_size, cell_start, cell_end, slice,
              forcing, u, rhs, lhsc, backsub_info, in_buffer, out_buffer, 
              #cv, rhon, rhos, rhoq, 
              cuf, q, ue, buf, fjac, njac, lhsa, lhsb, us, vs, ws, qs, rho_i, square
end

