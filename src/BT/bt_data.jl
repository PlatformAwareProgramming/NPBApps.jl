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

#const predecessor = zeros(Int64, 3)
#const successor = zeros(Int64,3)
#const grid_size = zeros(Int64, 3)

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
const t_rhsx = 14
const t_rhsy = 15
const t_rhsz = 16
const t_add = 17
const t_rdis1 = 18
const t_rdis2 = 19
const t_last = 19

const t_names = ["total", "i/o", "rhs", "xsolve", "ysolve", "zsolve",
                "bpack", "exch", "xcomm", "ycomm", "zcomm",
                " totcomp", " totcomm", "rhsx", "rhsy", "rhsz",  
                "add", "qbc_send", "qbc_recv", "last"]

t_recs = "total"

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

#---------------------------------------------------------------------
# allocate space dynamically for data arrays
#---------------------------------------------------------------------

 function alloc_field_space(z, grid_points, x_zones, y_zones)

      @info "grid_points=$grid_points"
      problem_size = maximum(div.(grid_points, [x_zones, y_zones, 1]) .+ 1)
      @info "grid_points=$grid_points, problem_size = $problem_size"

      MAX_CELL_DIM[z] = div(problem_size, maxcells)+1

#      IMAX[z] = MAX_CELL_DIM[z]
#      JMAX[z] = MAX_CELL_DIM[z]
#      KMAX[z] = MAX_CELL_DIM[z]

      IMAX[z] = div(grid_points[1], maxcells) + 1 #MAX_CELL_DIM
      JMAX[z] = div(grid_points[2], maxcells) + 1 #MAX_CELL_DIM
      KMAX[z] = div(grid_points[3], maxcells) + 1 #MAX_CELL_DIM


      BUF_SIZE[z] = MAX_CELL_DIM[z]*MAX_CELL_DIM[z]*(maxcells-1)*60+1
 #     @info "$clusterid/$node: **** BUF_SIZE[$z]=$(BUF_SIZE[z])  --- MAX_CELL_DIM[$z] = $(MAX_CELL_DIM[z]) --- maxcells=$maxcells"

      cell_coord[z] = zeros(Int64, 3, maxcells)
      cell_low[z] = zeros(Int64, 3, maxcells)
      cell_high[z] = zeros(Int64, 3, maxcells)
      cell_size[z] = zeros(Int64, 3, maxcells)
      cell_start[z] = zeros(Int64, 3, maxcells)
      cell_end[z] = zeros(Int64, 3, maxcells)
      slice[z] = zeros(Int64, 3, maxcells)

      predecessor[z] = zeros(Int64, 3)
      successor[z] = zeros(Int64,3)

      forcing0 = zeros(Float64, 5, IMAX[z], JMAX[z], KMAX[z], maxcells)
      forcing[z] = OffsetArray(forcing0, 1:5, 0:IMAX[z]-1, 0:JMAX[z]-1, 0:KMAX[z]-1, 1:maxcells)            

      u0 = zeros(Float64, 5, IMAX[z]+4, JMAX[z]+4, KMAX[z] + 4, maxcells)
      u[z] = OffsetArray(u0, 1:5, -2:IMAX[z]+1, -2:JMAX[z]+1, -2:KMAX[z]+1, 1:maxcells)

      rhs0 = zeros(Float64, 5, IMAX[z]+1, JMAX[z]+1, KMAX[z]+1, maxcells)
      rhs[z] = OffsetArray(rhs0, 1:5, -1:IMAX[z]-1, -1:JMAX[z]-1, -1:KMAX[z]-1, 1:maxcells)

      lhsc0 = zeros(Float64, 5, 5, IMAX[z]+1, JMAX[z]+1, KMAX[z]+1, maxcells)
      lhsc[z] = OffsetArray(lhsc0, 1:5, 1:5, -1:IMAX[z]-1, -1:JMAX[z]-1, -1:KMAX[z]-1, 1:maxcells)

      backsub_info0 = zeros(Float64, 5, MAX_CELL_DIM[z]+1, MAX_CELL_DIM[z]+1, maxcells)
      backsub_info[z] = OffsetArray(backsub_info0, 1:5, 0:MAX_CELL_DIM[z], 0:MAX_CELL_DIM[z], 1:maxcells)

      in_buffer[z] = Array{Float64}(undef, BUF_SIZE[z])
      out_buffer[z] = Array{Float64}(undef, BUF_SIZE[z])

      cv0 = zeros(Float64, MAX_CELL_DIM[z] + 4)
      cv[z] = OffsetArray(cv0, -2:MAX_CELL_DIM[z]+1)

      rhon0 = zeros(Float64, MAX_CELL_DIM[z] + 4)
      rhon[z] = OffsetArray(rhon0, -2:MAX_CELL_DIM[z]+1)

      rhos0 = zeros(Float64, MAX_CELL_DIM[z] + 4)
      rhos[z] = OffsetArray(rhos0, -2:MAX_CELL_DIM[z]+1)
      
      rhoq0 = zeros(Float64, MAX_CELL_DIM[z] + 4)
      rhoq[z] = OffsetArray(rhoq0, -2:MAX_CELL_DIM[z]+1)

      cuf0 = zeros(Float64, MAX_CELL_DIM[z] + 4)
      cuf[z] = OffsetArray(cuf0, -2:MAX_CELL_DIM[z]+1)

      q0 = zeros(Float64, MAX_CELL_DIM[z] + 4)
      q[z] = OffsetArray(q0, -2:MAX_CELL_DIM[z]+1)
      
      ue0 = zeros(Float64, MAX_CELL_DIM[z] + 4, 5)
      ue[z] = OffsetArray(ue0, -2:MAX_CELL_DIM[z]+1, 1:5)

      buf0 = zeros(Float64, MAX_CELL_DIM[z] + 4, 5)
      buf[z] = OffsetArray(buf0, -2:MAX_CELL_DIM[z]+1, 1:5)

      fjac0 = zeros(Float64, 5, 5, MAX_CELL_DIM[z] + 4)
      fjac[z] = OffsetArray(fjac0, 1:5, 1:5, -2:MAX_CELL_DIM[z]+1)

      njac0 = zeros(Float64, 5, 5, MAX_CELL_DIM[z] + 4)
      njac[z] = OffsetArray(njac0, 1:5, 1:5, -2:MAX_CELL_DIM[z]+1)

      lhsa0 = zeros(Float64, 5, 5, MAX_CELL_DIM[z] + 2)
      lhsa[z] = OffsetArray(lhsa0, 1:5, 1:5, -1:MAX_CELL_DIM[z])

      lhsb0 = zeros(Float64, 5, 5, MAX_CELL_DIM[z] + 2)
      lhsb[z] = OffsetArray(lhsb0, 1:5, 1:5, -1:MAX_CELL_DIM[z])            

      us0 = zeros(Float64, IMAX[z]+2, JMAX[z]+2, KMAX[z]+2, maxcells)
      us[z] = OffsetArray(us0, -1:IMAX[z], -1:JMAX[z], -1:KMAX[z], 1:maxcells)
      
      vs0 = zeros(Float64, IMAX[z]+2, JMAX[z]+2, KMAX[z]+2, maxcells)
      vs[z] = OffsetArray(vs0, -1:IMAX[z], -1:JMAX[z], -1:KMAX[z], 1:maxcells)

      ws0 = zeros(Float64, IMAX[z]+2, JMAX[z]+2, KMAX[z]+2, maxcells)
      ws[z] = OffsetArray(ws0, -1:IMAX[z], -1:JMAX[z], -1:KMAX[z], 1:maxcells)

      qs0 = zeros(Float64, IMAX[z]+2, JMAX[z]+2, KMAX[z]+2, maxcells)
      qs[z] = OffsetArray(qs0, -1:IMAX[z], -1:JMAX[z], -1:KMAX[z], 1:maxcells)
      
      rho_i0 = zeros(Float64, IMAX[z]+2, JMAX[z]+2, KMAX[z]+2, maxcells)
      rho_i[z] = OffsetArray(rho_i0, -1:IMAX[z], -1:JMAX[z], -1:KMAX[z], 1:maxcells)

      square0 = zeros(Float64, IMAX[z]+2, JMAX[z]+2, KMAX[z]+2, maxcells)
      square[z] = OffsetArray(square0, -1:IMAX[z], -1:JMAX[z], -1:KMAX[z], 1:maxcells)
      
      return nothing
end

function alloc_field_space_zones(proc_num_zones)

      global MAX_CELL_DIM = Array{Int64}(undef, proc_num_zones) 

      global IMAX = Array{Int64}(undef, proc_num_zones) 
      global JMAX = Array{Int64}(undef, proc_num_zones) 
      global KMAX = Array{Int64}(undef, proc_num_zones) 

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
      
      # field arrays

      global forcing = Array{OffsetArray{Float64, 5, Array{Float64, 5}}}(undef, proc_num_zones)
      global u = Array{OffsetArray{Float64, 5, Array{Float64, 5}}}(undef, proc_num_zones)
      global rhs = Array{OffsetArray{Float64, 5, Array{Float64, 5}}}(undef, proc_num_zones)
      global lhsc = Array{OffsetArray{Float64, 6, Array{Float64, 6}}}(undef, proc_num_zones)
      global backsub_info = Array{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, proc_num_zones)
      global in_buffer = Array{Array{Float64}}(undef, proc_num_zones)
      global out_buffer = Array{Array{Float64}}(undef, proc_num_zones)
      global cv = Array{OffsetArray{Float64, 1, Array{Float64, 1}}}(undef, proc_num_zones)
      global rhon = Array{OffsetArray{Float64, 1, Array{Float64, 1}}}(undef, proc_num_zones)
      global rhos = Array{OffsetArray{Float64, 1, Array{Float64, 1}}}(undef, proc_num_zones)
      global rhoq = Array{OffsetArray{Float64, 1, Array{Float64, 1}}}(undef, proc_num_zones)
      global cuf = Array{OffsetArray{Float64, 1, Array{Float64, 1}}}(undef, proc_num_zones)
      global q = Array{OffsetArray{Float64, 1, Array{Float64, 1}}}(undef, proc_num_zones)
      global ue = Array{OffsetArray{Float64, 2, Array{Float64, 2}}}(undef, proc_num_zones)
      global buf = Array{OffsetArray{Float64, 2, Array{Float64, 2}}}(undef, proc_num_zones)
      global fjac = Array{OffsetArray{Float64, 3, Array{Float64, 3}}}(undef, proc_num_zones)
      global njac = Array{OffsetArray{Float64, 3, Array{Float64, 3}}}(undef, proc_num_zones)
      global lhsa = Array{OffsetArray{Float64, 3, Array{Float64, 3}}}(undef, proc_num_zones)
      global lhsb = Array{OffsetArray{Float64, 3, Array{Float64, 3}}}(undef, proc_num_zones)
      global us = Array{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, proc_num_zones)
      global vs = Array{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, proc_num_zones)
      global ws = Array{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, proc_num_zones)
      global qs = Array{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, proc_num_zones)
      global ainv = Array{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, proc_num_zones)
      global rho_i = Array{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, proc_num_zones)
      global square = Array{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, proc_num_zones)
             
      return nothing

end

function print_u(z)

      for c = 1:ncells
            for kk = -1:KMAX[z]
               for jj = -1:JMAX[z]
                  for ii = -1:IMAX[z]
                        for m = 1:5
                            @info  "u[$z][$m, $ii, $jj, $kk, $c] = $(u[z][m, ii, jj, kk, c])"
                        end
                  end
               end
            end
      end

end