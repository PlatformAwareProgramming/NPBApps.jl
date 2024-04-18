#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---  lu_data module
#---------------------------------------------------------------------
#---------------------------------------------------------------------

const ipr_default = 1
const omega_default = 1.2e0
tolrsddef = Array{Float64}(undef, 5)
tolrsddef .= 1.0f-08

const c1 = 1.40e+00 
const c2 = 0.40e+00
const c3 = 1.00e-01 
const c4 = 1.00e+00
const c5 = 1.40e+00


const from_s = 1
const from_n = 2
const from_e = 3
const from_w = 4


#---------------------------------------------------------------------
#     Timer constants
#---------------------------------------------------------------------

const t_total = 1 
const t_rhs = 2 
const t_blts = 3 
const t_buts = 4 
const t_jacld = 5
const t_jacu = 6 
const t_exch = 7 
const t_lcomm = 8 
const t_ucomm = 9 
const t_rcomm = 10
const t_rdis1 = 11
const t_rdis2 = 12
const t_last = 12

const t_names = "total", "rhs", "blts", "buts", "#jacld", "#jacu",
               "exch", "lcomm", "ucomm", "rcomm", "qbc_send", "qbc_recv",
               " totcomp", " totcomm"

function alloc_zone_vectors(x_zones, y_zones)

      global nx = zeros(Int64, max_zones)
      global ny = zeros(Int64, max_zones)
      global nz = zeros(Int64, max_zones)
      global nxmax = zeros(Int64, max_zones)
      global nx1 = zeros(Int64, max_zones)

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

 function alloc_field_space(z, nx0, ny0, nz0, nx, ny, nz)

      rsdnm[z] = Array{Float64}(undef,5)
      errnm[z] = Array{Float64}(undef,5)

      u0 = zeros(Float64, 5, nx0 + 4, ny0 + 4, nz0)
      u[z] = OffsetArray(u0, 1:5, -1:nx0+2, -1:ny0+2, 1:nz0)

      rsd0 = zeros(Float64, 5, nx0 + 4, ny0 + 4, nz0)
      rsd[z] = OffsetArray(rsd0, 1:5, -1:nx0+2, -1:ny0+2, 1:nz0)

      frct0 = zeros(Float64, 5, nx0 + 4, ny0 + 4, nz0)
      frct[z] = OffsetArray(frct0, 1:5, -1:nx0+2, -1:ny0+2, 1:nz0)

      flux0 = zeros(Float64, 5, nx0 + 2, ny0 + 2, nz0)
      flux[z] = OffsetArray(flux0, 1:5, 0:nx0+1,  0:ny0+1, 1:nz0)

      a[z] = zeros(Float64, 5, 5, nx0)
      b[z] = zeros(Float64, 5, 5, nx0)
      c[z] = zeros(Float64, 5, 5, nx0)
      d[z] = zeros(Float64, 5, 5, nx0)

      phi10 = zeros(Float64, nx0 + 2, ny0 + 2)
      phi1[z] = OffsetArray(phi10, 0:nx0+1, 0:ny0+1)

      phi20 = zeros(Float64, nx0 + 2, ny0 + 2)
      phi2[z] = OffsetArray(phi20, 0:nx0+1, 0:ny0+1)

      buf[z] = zeros(Float64, 5, 2*nx0*ny0)
      buf1[z] = zeros(Float64, 5, 2*nx0*ny0)

      buf_exch_w_in[z] = Array{Float64}(undef, 5*nx*nz)
      buf_exch_e_in[z] = Array{Float64}(undef, 5*nx*nz)
      buf_exch_n_in[z] = Array{Float64}(undef, 5*ny*nz)
      buf_exch_s_in[z] = Array{Float64}(undef, 5*ny*nz)

      buf_exch_w_out[z] = Array{Float64}(undef, 5*nx*nz)
      buf_exch_e_out[z] = Array{Float64}(undef, 5*nx*nz)
      buf_exch_n_out[z] = Array{Float64}(undef, 5*ny*nz)
      buf_exch_s_out[z] = Array{Float64}(undef, 5*ny*nz)

      return nothing
      
end

function alloc_field_space_zones(proc_num_zones)

      global rsdnm = Array{Array{Float64}}(undef, proc_num_zones)
      global errnm = Array{Array{Float64}}(undef, proc_num_zones)

      global u = Array{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, proc_num_zones)
      global rsd = Array{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, proc_num_zones)
      global frct = Array{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, proc_num_zones)
      global flux = Array{OffsetArray{Float64, 4, Array{Float64, 4}}}(undef, proc_num_zones)

      global a = Array{Array{Float64, 3}}(undef, proc_num_zones)
      global b = Array{Array{Float64, 3}}(undef, proc_num_zones) 
      global c = Array{Array{Float64, 3}}(undef, proc_num_zones) 
      global d = Array{Array{Float64, 3}}(undef, proc_num_zones) 

      global phi1 = Array{OffsetArray{Float64, 2, Array{Float64, 2}}}(undef, proc_num_zones)
      global phi2 = Array{OffsetArray{Float64, 2, Array{Float64, 2}}}(undef, proc_num_zones)

      global buf = Array{Array{Float64, 2}}(undef, proc_num_zones) 
      global buf1 = Array{Array{Float64, 2}}(undef, proc_num_zones) 

      global buf_exch_w_in = Array{Array{Float64}}(undef, proc_num_zones) 
      global buf_exch_e_in = Array{Array{Float64}}(undef, proc_num_zones) 
      global buf_exch_n_in = Array{Array{Float64}}(undef, proc_num_zones) 
      global buf_exch_s_in = Array{Array{Float64}}(undef, proc_num_zones) 

      global buf_exch_w_out = Array{Array{Float64}}(undef, proc_num_zones) 
      global buf_exch_e_out = Array{Array{Float64}}(undef, proc_num_zones) 
      global buf_exch_n_out = Array{Array{Float64}}(undef, proc_num_zones) 
      global buf_exch_s_out = Array{Array{Float64}}(undef, proc_num_zones) 

      return nothing
            
end