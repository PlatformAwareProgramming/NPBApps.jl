#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#  sp_data module
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------

#module sp_data

#---------------------------------------------------------------------
# The following include file is generated automatically by the
# "setparams" utility. It defines 
#      maxcells:      the square root of the maximum number of processors
#      problem_size:  12, 64, 102, 162 (for class S, A, B, C)
#      dt_default:    default time step for this problem size if no
#                     config file
#      niter_default: default number of iterations for this problem size
#---------------------------------------------------------------------

grid_points = zeros(Integer, 3)

const EAST = 2000 
const WEST = 3000      
const NORTH = 4000 
const SOUTH = 5000
const BOTTOM = 6000 
const TOP = 7000

predecessor = zeros(Int64, 3)
sucessor = zeros(Int64,3)
grid_size = zeros(Int64, 3)


#---------------------------------------------------------------------
#     Timer constants
#---------------------------------------------------------------------

const t_total = 1 
const t_rhs = 2 
const t_xsolve = 3 
const t_ysolve = 4
const t_zsolve = 5 
const t_bpack = 6 
const t_exch = 7 
const t_xcomm = 8
const t_ycomm = 9
const t_zcomm = 10 
const t_last = 10


#---------------------------------------------------------------------
#---------------------------------------------------------------------

function alloc_space()

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# allocate space dynamically for data arrays
#---------------------------------------------------------------------

      global MAX_CELL_DIM = div(problem_size,maxcells)+1

      global IMAX = MAX_CELL_DIM
      global JMAX = MAX_CELL_DIM
      global KMAX = MAX_CELL_DIM

      global IMAXP = div(IMAX,2)*2+1
      global JMAXP = div(JMAX,2)*2+1

#---------------------------------------------------------------------
# +1 at end to avoid zero length arrays for 1 node
#---------------------------------------------------------------------
      global BUF_SIZE = MAX_CELL_DIM*MAX_CELL_DIM*(maxcells-1)*60*2+1

      global cell_coord = zeros(Int64, 3, maxcells)
      global cell_low = zeros(Int64, 3, maxcells)
      global cell_high = zeros(Int64, 3, maxcells)
      global cell_size = zeros(Int64, 3, maxcells)
      global cell_start = zeros(Int64, 3, maxcells)
      global cell_end = zeros(Int64, 3, maxcells)
      global slice = zeros(Int64, 3, maxcells)

      u0 = zeros(Float64, IMAXP+4, JMAXP+4, KMAX + 4, 5, maxcells)
      global u = OffsetArray(u0, -2:IMAXP+1, -2:JMAXP+1, -2:KMAX+1, 1:5, 1:maxcells)

      us0 = zeros(Float64, IMAX+2, JMAX+2, KMAX+2, maxcells)
      global us = OffsetArray(us0, -1:IMAX, -1:JMAX, -1:KMAX, 1:maxcells)
      
      vs0 = zeros(Float64, IMAX+2, JMAX+2, KMAX+2, maxcells)
      global vs = OffsetArray(vs0, -1:IMAX, -1:JMAX, -1:KMAX, 1:maxcells)

      ws0 = zeros(Float64, IMAX+2, JMAX+2, KMAX+2, maxcells)
      global ws = OffsetArray(ws0, -1:IMAX, -1:JMAX, -1:KMAX, 1:maxcells)

      qs0 = zeros(Float64, IMAX+2, JMAX+2, KMAX+2, maxcells)
      global qs = OffsetArray(qs0, -1:IMAX, -1:JMAX, -1:KMAX, 1:maxcells)
     
      ainv0 = zeros(Float64, IMAX+2, JMAX+2, KMAX+2, maxcells)
      global ainv = OffsetArray(ainv0, -1:IMAX, -1:JMAX, -1:KMAX, 1:maxcells)
     
      rho_i0 = zeros(Float64, IMAX+2, JMAX+2, KMAX+2, maxcells)
      global rho_i = OffsetArray(rho_i0, -1:IMAX, -1:JMAX, -1:KMAX, 1:maxcells)

      speed0 = zeros(Float64, IMAX+2, JMAX+2, KMAX+2, maxcells)
      global speed = OffsetArray(speed0, -1:IMAX, -1:JMAX, -1:KMAX, 1:maxcells)
      
      square0 = zeros(Float64, IMAX+2, JMAX+2, KMAX+2, maxcells)
      global square = OffsetArray(square0, -1:IMAX, -1:JMAX, -1:KMAX, 1:maxcells)

      rhs0 = zeros(Float64, IMAXP, JMAXP, KMAX, 5, maxcells)
      global rhs = OffsetArray(rhs0, 0:IMAXP-1, 0:JMAXP-1, 0:KMAX-1, 1:5, 1:maxcells)

      forcing0 = zeros(Float64, IMAXP, JMAXP, KMAX, 5, maxcells)
      global forcing = OffsetArray(forcing0, 0:IMAXP-1, 0:JMAXP-1, 0:KMAX-1, 1:5, 1:maxcells)
      
      lhs0 = zeros(Float64, IMAXP, JMAXP, KMAX, 15, maxcells)
      global lhs = OffsetArray(lhs0, 0:IMAXP-1, 0:JMAXP-1, 0:KMAX-1, 1:15, 1:maxcells)

      global in_buffer = Array{Float64}(undef, BUF_SIZE)
      global out_buffer = Array{Float64}(undef, BUF_SIZE)

      cv0 = zeros(Float64, MAX_CELL_DIM + 4)
      global cv = OffsetArray(cv0, -2:MAX_CELL_DIM+1)

      rhon0 = zeros(Float64, MAX_CELL_DIM + 4)
      global rhon = OffsetArray(rhon0, -2:MAX_CELL_DIM+1)

      rhos0 = zeros(Float64, MAX_CELL_DIM + 4)
      global rhos = OffsetArray(rhos0, -2:MAX_CELL_DIM+1)
      
      rhoq0 = zeros(Float64, MAX_CELL_DIM + 4)
      global rhoq = OffsetArray(rhoq0, -2:MAX_CELL_DIM+1)

      cuf0 = zeros(Float64, MAX_CELL_DIM + 4)
      global cuf = OffsetArray(cuf0, -2:MAX_CELL_DIM+1)

      q0 = zeros(Float64, MAX_CELL_DIM + 4)
      global q = OffsetArray(q0, -2:MAX_CELL_DIM+1)
      
      ue0 = zeros(Float64, MAX_CELL_DIM + 4, 5)
      global ue = OffsetArray(ue0, -2:MAX_CELL_DIM+1, 1:5)

      buf0 = zeros(Float64, MAX_CELL_DIM + 4, 5)
      global buf = OffsetArray(buf0, -2:MAX_CELL_DIM+1, 1:5)
               
      return nothing
end

