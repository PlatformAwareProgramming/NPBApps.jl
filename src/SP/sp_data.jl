#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#  sp_data module
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

const EAST = 2000 
const WEST = 3000      
const NORTH = 4000 
const SOUTH = 5000
const BOTTOM = 6000 
const TOP = 7000

#---------------------------------------------------------------------
#     Timer constants
#---------------------------------------------------------------------

const t_total = 1
const t_rhsx = 2
const t_rhsy = 3
const t_rhsz = 4
const t_rhs = 5
const t_xsolve = 6
const t_ysolve = 7
const t_zsolve = 8
const t_txinvr = 9
const t_pinvr = 10
const t_ninvr = 11
const t_tzetar = 12
const t_add = 13
const t_rdis1 = 14
const t_rdis2 = 15
const t_bpack = 16 
const t_exch = 17 
const t_xcomm = 18
const t_ycomm = 19
const t_zcomm = 20 
const t_last = 20



function alloc_zone_vectors(x_zones, y_zones)

      global nx = zeros(Int64, max_zones)
      global ny = zeros(Int64, max_zones)
      global nz = zeros(Int64, max_zones)
      global nxmax = zeros(Int64, max_zones)

      global x_start = zeros(Int64, x_zones)   
      global y_start = zeros(Int64, y_zones)   
      global x_end = zeros(Int64, x_zones)     
      global y_end = zeros(Int64, y_zones)     
      global x_size = zeros(Int64, x_zones)    
      global y_size = zeros(Int64, y_zones)    

      global iz_west = zeros(Int64, max_zones)    
      global iz_east = zeros(Int64, max_zones)    
      global iz_south = zeros(Int64, max_zones)   
      global iz_north = zeros(Int64, max_zones)   

 
end


total_size = Ref{Int64}(0)

#---------------------------------------------------------------------
# allocate space dynamically for data arrays
#---------------------------------------------------------------------

function alloc_field_space(z, grid_points, x_zones, y_zones)

      problem_size = maximum(grid_points)

      MAX_CELL_DIM[z] = div(problem_size, maxcells) + 1

      IMAX[z] = div(grid_points[1], maxcells) + 1 #MAX_CELL_DIM
      JMAX[z] = div(grid_points[2], maxcells) + 1 #MAX_CELL_DIM
      KMAX[z] = div(grid_points[3], maxcells) + 1 #MAX_CELL_DIM

      IMAXP[z] = div(IMAX[z],2)*2+1
      JMAXP[z] = div(JMAX[z],2)*2+1

#---------------------------------------------------------------------
# +1 at end to avoid zero length arrays for 1 node
#---------------------------------------------------------------------
      BUF_SIZE[z] = MAX_CELL_DIM[z]*MAX_CELL_DIM[z]*(maxcells)*60*2+1

      cell_coord[z] = zeros(Int64, 3, maxcells)
      cell_low[z] = zeros(Int64, 3, maxcells)
      cell_high[z] = zeros(Int64, 3, maxcells)
      cell_size[z] = zeros(Int64, 3, maxcells)
      cell_start[z] = zeros(Int64, 3, maxcells)
      cell_end[z] = zeros(Int64, 3, maxcells)
      slice[z] = zeros(Int64, 3, maxcells)
      predecessor[z] = zeros(Int64, 3)
      successor[z] = zeros(Int64,3)
#      grid_size[z] = zeros(Int64, 3)
      
      # field arrays

      a = 0

      u0 = zeros(Float64, IMAXP[z]+4, JMAXP[z]+4, KMAX[z] + 4, 5, maxcells)
      u[z] = OffsetArray(u0, -2:IMAXP[z]+1, -2:JMAXP[z]+1, -2:KMAX[z]+1, 1:5, 1:maxcells)

      total_size[] += sizeof(u[z])

      us0 = zeros(Float64, IMAX[z]+2, JMAX[z]+2, KMAX[z]+2, maxcells)
      us[z] = OffsetArray(us0, -1:IMAX[z], -1:JMAX[z], -1:KMAX[z], 1:maxcells)
      
      total_size[] += sizeof(us[z])

      vs0 = zeros(Float64, IMAX[z]+2, JMAX[z]+2, KMAX[z]+2, maxcells)
      vs[z] = OffsetArray(vs0, -1:IMAX[z], -1:JMAX[z], -1:KMAX[z], 1:maxcells)

      total_size[] += sizeof(vs[z])

      ws0 = zeros(Float64, IMAX[z]+2, JMAX[z]+2, KMAX[z]+2, maxcells)
      ws[z] = OffsetArray(ws0, -1:IMAX[z], -1:JMAX[z], -1:KMAX[z], 1:maxcells)

      total_size[] += sizeof(ws[z])

      qs0 = zeros(Float64, IMAX[z]+2, JMAX[z]+2, KMAX[z]+2, maxcells)
      qs[z] = OffsetArray(qs0, -1:IMAX[z], -1:JMAX[z], -1:KMAX[z], 1:maxcells)
     
      total_size[] += sizeof(qs[z])

      ainv0 = zeros(Float64, IMAX[z]+2, JMAX[z]+2, KMAX[z]+2, maxcells)
      ainv[z] = OffsetArray(ainv0, -1:IMAX[z], -1:JMAX[z], -1:KMAX[z], 1:maxcells)
     
      total_size[] += sizeof(ainv[z])

      rho_i0 = zeros(Float64, IMAX[z]+2, JMAX[z]+2, KMAX[z]+2, maxcells)
      rho_i[z] = OffsetArray(rho_i0, -1:IMAX[z], -1:JMAX[z], -1:KMAX[z], 1:maxcells)

      total_size[] += sizeof(rho_i[z])

      speed0 = zeros(Float64, IMAX[z]+2, JMAX[z]+2, KMAX[z]+2, maxcells)
      speed[z] = OffsetArray(speed0, -1:IMAX[z], -1:JMAX[z], -1:KMAX[z], 1:maxcells)
      
      total_size[] += sizeof(speed[z])

      square0 = zeros(Float64, IMAX[z]+2, JMAX[z]+2, KMAX[z]+2, maxcells)
      square[z] = OffsetArray(square0, -1:IMAX[z], -1:JMAX[z], -1:KMAX[z], 1:maxcells)

      total_size[] += sizeof(square[z])

      rhs0 = zeros(Float64, IMAXP[z], JMAXP[z], KMAX[z], 5, maxcells)
      rhs[z] = OffsetArray(rhs0, 0:IMAXP[z]-1, 0:JMAXP[z]-1, 0:KMAX[z]-1, 1:5, 1:maxcells)

      total_size[] += sizeof(rhs[z])

      forcing0 = zeros(Float64, IMAXP[z], JMAXP[z], KMAX[z], 5, maxcells)
      forcing[z] = OffsetArray(forcing0, 0:IMAXP[z]-1, 0:JMAXP[z]-1, 0:KMAX[z]-1, 1:5, 1:maxcells)
      
      total_size[] += sizeof(forcing[z])

      lhs0 = zeros(Float64, IMAXP[z], JMAXP[z], KMAX[z], 15, maxcells)
      lhs[z] = OffsetArray(lhs0, 0:IMAXP[z]-1, 0:JMAXP[z]-1, 0:KMAX[z]-1, 1:15, 1:maxcells)

      total_size[] += sizeof(lhs[z])

      in_buffer[z] = Array{Float64}(undef, BUF_SIZE[z])
      out_buffer[z] = Array{Float64}(undef, BUF_SIZE[z])

      total_size[] += sizeof(in_buffer[z])
      total_size[] += sizeof(out_buffer[z])

      cv0 = zeros(Float64, MAX_CELL_DIM[z] + 4)
      cv[z] = OffsetArray(cv0, -2:MAX_CELL_DIM[z]+1)

      total_size[] += sizeof(cv[z])

      rhon0 = zeros(Float64, MAX_CELL_DIM[z] + 4)
      rhon[z] = OffsetArray(rhon0, -2:MAX_CELL_DIM[z]+1)

      total_size[] += sizeof(rhon[z])
      
      rhos0 = zeros(Float64, MAX_CELL_DIM[z] + 4)
      rhos[z] = OffsetArray(rhos0, -2:MAX_CELL_DIM[z]+1)
      
      total_size[] += sizeof(rhos[z])
      
      rhoq0 = zeros(Float64, MAX_CELL_DIM[z] + 4)
      rhoq[z] = OffsetArray(rhoq0, -2:MAX_CELL_DIM[z]+1)

      total_size[] += sizeof(rhoq[z])
      
      cuf0 = zeros(Float64, MAX_CELL_DIM[z] + 4)
      cuf[z] = OffsetArray(cuf0, -2:MAX_CELL_DIM[z]+1)

      total_size[] += sizeof(cuf[z])
      
      q0 = zeros(Float64, MAX_CELL_DIM[z] + 4)
      q[z] = OffsetArray(q0, -2:MAX_CELL_DIM[z]+1)
      
      ue0 = zeros(Float64, MAX_CELL_DIM[z] + 4, 5)
      ue[z] = OffsetArray(ue0, -2:MAX_CELL_DIM[z]+1, 1:5)

      total_size[] += sizeof(ue[z])
      
      buf0 = zeros(Float64, MAX_CELL_DIM[z] + 4, 5)
      buf[z] = OffsetArray(buf0, -2:MAX_CELL_DIM[z]+1, 1:5)               
      
      total_size[] += sizeof(buf[z])
      
      return nothing
end


function alloc_field_space_zones(proc_num_zones)

      global MAX_CELL_DIM = Array{Int64}(undef, proc_num_zones) 

      global IMAX = Array{Int64}(undef, proc_num_zones) 
      global JMAX = Array{Int64}(undef, proc_num_zones) 
      global KMAX = Array{Int64}(undef, proc_num_zones) 

      global IMAXP = Array{Int64}(undef, proc_num_zones)
      global JMAXP = Array{Int64}(undef, proc_num_zones) 

      global BUF_SIZE = Array{Int64}(undef, proc_num_zones) 

#---------------------------------------------------------------------
# +1 at end to avoid zero length arrays for 1 node
#---------------------------------------------------------------------

      global cell_coord = Array{Array{Int64,2}}(undef, proc_num_zones) 
      global cell_low = Array{Array{Int64,2}}(undef, proc_num_zones) 
      global cell_high = Array{Array{Int64,2}}(undef, proc_num_zones) 
      global cell_size = Array{Array{Int64,2}}(undef, proc_num_zones) 
      global cell_start = Array{Array{Int64,2}}(undef, proc_num_zones) 
      global cell_end = Array{Array{Int64,2}}(undef, proc_num_zones) 
      global slice = Array{Array{Int64,2}}(undef, proc_num_zones) 
      global predecessor = Array{Array{Int64}}(undef, proc_num_zones) 
      global successor = Array{Array{Int64}}(undef, proc_num_zones) 
#      global grid_size = Array{Array{Int64}}(undef, proc_num_zones) 
      
      # field arrays

      global u = Array{OffsetArray{Float64, 5, Array{Float64, 5}}}(undef, proc_num_zones)
      global us = Array{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, proc_num_zones)
      global vs = Array{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, proc_num_zones)
      global ws = Array{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, proc_num_zones)
      global qs = Array{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, proc_num_zones)
      global ainv = Array{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, proc_num_zones)
      global rho_i = Array{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, proc_num_zones)
      global speed = Array{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, proc_num_zones)
      global square = Array{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, proc_num_zones)
      global rhs = Array{OffsetArray{Float64, 5, Array{Float64, 5}}}(undef, proc_num_zones)
      global forcing = Array{OffsetArray{Float64, 5, Array{Float64, 5}}}(undef, proc_num_zones)
      global lhs = Array{OffsetArray{Float64, 5, Array{Float64, 5}}}(undef, proc_num_zones)
      global cv = Array{OffsetArray{Float64, 1, Array{Float64, 1}}}(undef, proc_num_zones)
      global rhon = Array{OffsetArray{Float64, 1, Array{Float64, 1}}}(undef, proc_num_zones)
      global rhos = Array{OffsetArray{Float64, 1, Array{Float64, 1}}}(undef, proc_num_zones)
      global rhoq = Array{OffsetArray{Float64, 1, Array{Float64, 1}}}(undef, proc_num_zones)
      global cuf = Array{OffsetArray{Float64, 1, Array{Float64, 1}}}(undef, proc_num_zones)
      global q = Array{OffsetArray{Float64, 1, Array{Float64, 1}}}(undef, proc_num_zones)
      global ue = Array{OffsetArray{Float64, 2, Array{Float64, 2}}}(undef, proc_num_zones)
      global buf = Array{OffsetArray{Float64, 2, Array{Float64, 2}}}(undef, proc_num_zones)
               
      global in_buffer = Array{Array{Float64}}(undef, proc_num_zones)
      global out_buffer = Array{Array{Float64}}(undef, proc_num_zones)

      return nothing
end