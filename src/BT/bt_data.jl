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

 function alloc_field_space(z, grid_points, problem_size)

      total_size = Ref{Int64}(0)

#      MAX_CELL_DIM[z] = div(problem_size, maxcells)+1    #64
      @warn "------>>>>> zone=$z $problem_size $maxcells $(MAX_CELL_DIM[z])"

#      IMAX[z] = MAX_CELL_DIM[z]
#      JMAX[z] = MAX_CELL_DIM[z]
#      KMAX[z] = MAX_CELL_DIM[z]

      IMAX[z] = div(grid_points[1], maxcells) + 1 #MAX_CELL_DIM
      JMAX[z] = div(grid_points[2], maxcells) + 1 #MAX_CELL_DIM
      KMAX[z] = div(grid_points[3], maxcells) + 1 #MAX_CELL_DIM

      MAX_CELL_DIM[z] = max(IMAX[z], JMAX[z], KMAX[z])

      BUF_SIZE[z] = MAX_CELL_DIM[z]*MAX_CELL_DIM[z]*(maxcells)*60+1

      cell_coord[z] = zeros(Int64, 3, maxcells)
      cell_low[z] = zeros(Int64, 3, maxcells)
      cell_high[z] = zeros(Int64, 3, maxcells)
      cell_size[z] = zeros(Int64, 3, maxcells)
      cell_start[z] = zeros(Int64, 3, maxcells)
      cell_end[z] = zeros(Int64, 3, maxcells)
      slice[z] = zeros(Int64, 3, maxcells)

      predecessor[z] = zeros(Int64, 3)
      successor[z] = zeros(Int64,3)

      forcing0 = zeros(FloatType, 5, IMAX[z], JMAX[z], KMAX[z], maxcells)
      forcing[z] = OffsetArray(forcing0, 1:5, 0:IMAX[z]-1, 0:JMAX[z]-1, 0:KMAX[z]-1, 1:maxcells)    # 20,971,520

      total_size[] += sizeof(forcing[z])

      u0 = zeros(FloatType, 5, IMAX[z]+4, JMAX[z]+4, KMAX[z] + 4, maxcells)
      u[z] = OffsetArray(u0, 1:5, -2:IMAX[z]+1, -2:JMAX[z]+1, -2:KMAX[z]+1, 1:maxcells)   # 25,154,560

      total_size[] += sizeof(u[z])

      rhs0 = zeros(FloatType, 5, IMAX[z]+1, JMAX[z]+1, KMAX[z]+1, maxcells)
      rhs[z] = OffsetArray(rhs0, 1:5, -1:IMAX[z]-1, -1:JMAX[z]-1, -1:KMAX[z]-1, 1:maxcells) # 21.970.000

      total_size[] += sizeof(rhs[z])

      lhsc0 = zeros(FloatType, 5, 5, IMAX[z]+1, JMAX[z]+1, KMAX[z]+1, maxcells) # 21.970.000
      lhsc[z] = OffsetArray(lhsc0, 1:5, 1:5, -1:IMAX[z]-1, -1:JMAX[z]-1, -1:KMAX[z]-1, 1:maxcells)

      total_size[] += sizeof(lhsc[z])

      backsub_info0 = zeros(FloatType, 5, MAX_CELL_DIM[z]+1, MAX_CELL_DIM[z]+1, maxcells)
      backsub_info[z] = OffsetArray(backsub_info0, 1:5, 0:MAX_CELL_DIM[z], 0:MAX_CELL_DIM[z], 1:maxcells)

      total_size[] += sizeof(backsub_info[z])


      cv0 = zeros(FloatType, MAX_CELL_DIM[z] + 4)
      cv[z] = OffsetArray(cv0, -2:MAX_CELL_DIM[z]+1)

      total_size[] += sizeof(cv[z])

      rhon0 = zeros(FloatType, MAX_CELL_DIM[z] + 4)
      rhon[z] = OffsetArray(rhon0, -2:MAX_CELL_DIM[z]+1)

      total_size[] += sizeof(rhon[z])

      rhos0 = zeros(FloatType, MAX_CELL_DIM[z] + 4)
      rhos[z] = OffsetArray(rhos0, -2:MAX_CELL_DIM[z]+1)
      
      total_size[] += sizeof(rhos[z])

      rhoq0 = zeros(FloatType, MAX_CELL_DIM[z] + 4)
      rhoq[z] = OffsetArray(rhoq0, -2:MAX_CELL_DIM[z]+1)

      total_size[] += sizeof(rhoq[z])

      cuf0 = zeros(FloatType, MAX_CELL_DIM[z] + 4)
      cuf[z] = OffsetArray(cuf0, -2:MAX_CELL_DIM[z]+1)

      total_size[] += sizeof(cuf[z])

      q0 = zeros(FloatType, MAX_CELL_DIM[z] + 4)
      q[z] = OffsetArray(q0, -2:MAX_CELL_DIM[z]+1)
      
      total_size[] += sizeof(q[z])

      ue0 = zeros(FloatType, MAX_CELL_DIM[z] + 4, 5)
      ue[z] = OffsetArray(ue0, -2:MAX_CELL_DIM[z]+1, 1:5)

      total_size[] += sizeof(ue[z])

      buf0 = zeros(FloatType, MAX_CELL_DIM[z] + 4, 5)
      buf[z] = OffsetArray(buf0, -2:MAX_CELL_DIM[z]+1, 1:5)

      total_size[] += sizeof(buf[z])

      fjac0 = zeros(FloatType, 5, 5, MAX_CELL_DIM[z] + 4)
      fjac[z] = OffsetArray(fjac0, 1:5, 1:5, -2:MAX_CELL_DIM[z]+1)

      total_size[] += sizeof(fjac[z])

      njac0 = zeros(FloatType, 5, 5, MAX_CELL_DIM[z] + 4)
      njac[z] = OffsetArray(njac0, 1:5, 1:5, -2:MAX_CELL_DIM[z]+1)

      total_size[] += sizeof(njac[z])


   #   if z == 1

            in_buffer[z] = Array{FloatType}(undef, BUF_SIZE[z])
            out_buffer[z] = Array{FloatType}(undef, BUF_SIZE[z])
       
            total_size[] += sizeof(in_buffer[z])
            total_size[] += sizeof(out_buffer[z])
      
   #   else
   #         in_buffer[z] = in_buffer[1]
   #         out_buffer[z] = out_buffer[1]
   #   end

      lhsa0 = zeros(FloatType, 5, 5, MAX_CELL_DIM[z] + 2)
      lhsa[z] = OffsetArray(lhsa0, 1:5, 1:5, -1:MAX_CELL_DIM[z])

      total_size[] += sizeof(lhsa[z])

      lhsb0 = zeros(FloatType, 5, 5, MAX_CELL_DIM[z] + 2)
      lhsb[z] = OffsetArray(lhsb0, 1:5, 1:5, -1:MAX_CELL_DIM[z])            

      total_size[] += sizeof(lhsb[z])

      us0 = zeros(FloatType, IMAX[z]+2, JMAX[z]+2, KMAX[z]+2, maxcells)
      us[z] = OffsetArray(us0, -1:IMAX[z], -1:JMAX[z], -1:KMAX[z], 1:maxcells)      # 4,599,936
      total_size[] += sizeof(us[z])

      vs0 = zeros(FloatType, IMAX[z]+2, JMAX[z]+2, KMAX[z]+2, maxcells)
      vs[z] = OffsetArray(vs0, -1:IMAX[z], -1:JMAX[z], -1:KMAX[z], 1:maxcells)      # 4,599,936

      total_size[] += sizeof(vs[z])

      ws0 = zeros(FloatType, IMAX[z]+2, JMAX[z]+2, KMAX[z]+2, maxcells)
      ws[z] = OffsetArray(ws0, -1:IMAX[z], -1:JMAX[z], -1:KMAX[z], 1:maxcells)      # 4,599,936

      total_size[] += sizeof(ws[z])

      qs0 = zeros(FloatType, IMAX[z]+2, JMAX[z]+2, KMAX[z]+2, maxcells)
      qs[z] = OffsetArray(qs0, -1:IMAX[z], -1:JMAX[z], -1:KMAX[z], 1:maxcells)      # 4,599,936
      
      total_size[] += sizeof(qs[z])

      rho_i0 = zeros(FloatType, IMAX[z]+2, JMAX[z]+2, KMAX[z]+2, maxcells)
      rho_i[z] = OffsetArray(rho_i0, -1:IMAX[z], -1:JMAX[z], -1:KMAX[z], 1:maxcells)      # 4,599,936

      total_size[] += sizeof(rho_i[z])

      square0 = zeros(FloatType, IMAX[z]+2, JMAX[z]+2, KMAX[z]+2, maxcells)
      square[z] = OffsetArray(square0, -1:IMAX[z], -1:JMAX[z], -1:KMAX[z], 1:maxcells)      # 4,599,936
      
      total_size[] += sizeof(square[z])

      return total_size[]
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

      global forcing = Array{OffsetArray{FloatType, 5, Array{FloatType, 5}}}(undef, proc_num_zones)
      global u = Array{OffsetArray{FloatType, 5, Array{FloatType, 5}}}(undef, proc_num_zones)
      global rhs = Array{OffsetArray{FloatType, 5, Array{FloatType, 5}}}(undef, proc_num_zones)
      global lhsc = Array{OffsetArray{FloatType, 6, Array{FloatType, 6}}}(undef, proc_num_zones)
      global backsub_info = Array{OffsetArray{FloatType, 4, Array{FloatType, 4}}}(undef, proc_num_zones)
      global in_buffer = Array{Array{FloatType}}(undef, proc_num_zones)
      global out_buffer = Array{Array{FloatType}}(undef, proc_num_zones)
      global cv = Array{OffsetArray{FloatType, 1, Array{FloatType, 1}}}(undef, proc_num_zones)
      global rhon = Array{OffsetArray{FloatType, 1, Array{FloatType, 1}}}(undef, proc_num_zones)
      global rhos = Array{OffsetArray{FloatType, 1, Array{FloatType, 1}}}(undef, proc_num_zones)
      global rhoq = Array{OffsetArray{FloatType, 1, Array{FloatType, 1}}}(undef, proc_num_zones)
      global cuf = Array{OffsetArray{FloatType, 1, Array{FloatType, 1}}}(undef, proc_num_zones)
      global q = Array{OffsetArray{FloatType, 1, Array{FloatType, 1}}}(undef, proc_num_zones)
      global ue = Array{OffsetArray{FloatType, 2, Array{FloatType, 2}}}(undef, proc_num_zones)
      global buf = Array{OffsetArray{FloatType, 2, Array{FloatType, 2}}}(undef, proc_num_zones)
      global fjac = Array{OffsetArray{FloatType, 3, Array{FloatType, 3}}}(undef, proc_num_zones)
      global njac = Array{OffsetArray{FloatType, 3, Array{FloatType, 3}}}(undef, proc_num_zones)
      global lhsa = Array{OffsetArray{FloatType, 3, Array{FloatType, 3}}}(undef, proc_num_zones)
      global lhsb = Array{OffsetArray{FloatType, 3, Array{FloatType, 3}}}(undef, proc_num_zones)
      global us = Array{OffsetArray{FloatType, 4, Array{FloatType, 4}}}(undef, proc_num_zones)
      global vs = Array{OffsetArray{FloatType, 4, Array{FloatType, 4}}}(undef, proc_num_zones)
      global ws = Array{OffsetArray{FloatType, 4, Array{FloatType, 4}}}(undef, proc_num_zones)
      global qs = Array{OffsetArray{FloatType, 4, Array{FloatType, 4}}}(undef, proc_num_zones)
      global ainv = Array{OffsetArray{FloatType, 4, Array{FloatType, 4}}}(undef, proc_num_zones)
      global rho_i = Array{OffsetArray{FloatType, 4, Array{FloatType, 4}}}(undef, proc_num_zones)
      global square = Array{OffsetArray{FloatType, 4, Array{FloatType, 4}}}(undef, proc_num_zones)
             
      return nothing

end

function print_u(z)

      for c = 1:ncells
            for kk = -1:KMAX[z]
               for jj = -1:JMAX[z]
                  for ii = -1:IMAX[z]
                        for m = 1:5
#                            @info  "u[$z][$m, $ii, $jj, $kk, $c] = $(u[z][m, ii, jj, kk, c])"
                             @info  "$ii $jj $kk $m $c $(u[z][m, ii, jj, kk, c])"
                        end
                  end
               end
            end
      end

end