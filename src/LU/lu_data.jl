#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---  lu_data module
#---------------------------------------------------------------------
#---------------------------------------------------------------------

const ipr_default = 1
const omega_default = 1.2e0
tolrsddef = Array{FloatType}(undef, 5)
tolrsddef .= 1.0f-08

const c1 = 1.40e+00 
const c2 = 0.40e+00
const c3 = 1.00e-01 
const c4 = 1.00e+00
const c5 = 1.40e+00


const from_s = 1000
const from_n = 2000
const from_e = 3000
const from_w = 4000


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

 function alloc_field_space(z, nx, ny, nz, problem_size)

      rsdnm[z] = Array{FloatType}(undef,5)
      errnm[z] = Array{FloatType}(undef,5)

      u0 = zeros(FloatType, 5, nx + 4, ny + 4, nz)
      u[z] = OffsetArray(u0, 1:5, -1:nx+2, -1:ny+2, 1:nz)

      rsd0 = zeros(FloatType, 5, nx + 4, ny + 4, nz)
      rsd[z] = OffsetArray(rsd0, 1:5, -1:nx+2, -1:ny+2, 1:nz)

      frct0 = zeros(FloatType, 5, nx + 4, ny + 4, nz)
      frct[z] = OffsetArray(frct0, 1:5, -1:nx+2, -1:ny+2, 1:nz)

      flux0 = zeros(FloatType, 5, nx + 2, ny + 2, nz)
      flux[z] = OffsetArray(flux0, 1:5, 0:nx+1,  0:ny+1, 1:nz)

      a[z] = zeros(FloatType, 5, 5, nx)
      b[z] = zeros(FloatType, 5, 5, nx)
      c[z] = zeros(FloatType, 5, 5, nx)
      d[z] = zeros(FloatType, 5, 5, nx)

      nmax = max(nx, ny, nz)

      phi10 = zeros(FloatType, nmax + 2, nmax + 2)
      phi1[z] = OffsetArray(phi10, 0:nmax+1, 0:nmax+1)

      phi20 = zeros(FloatType, nmax + 2, nmax + 2)
      phi2[z] = OffsetArray(phi20, 0:nmax+1, 0:nmax+1)

      buf[z] = zeros(FloatType, 5, 2*max(nx,ny)*nz)
      buf1[z] = zeros(FloatType, 5, 2*max(nx,ny)*nz)

      buf_exch_w_in[z] = Array{FloatType}(undef, 5*(ny#=*-2*=#)*(nz-2))
      buf_exch_e_in[z] = Array{FloatType}(undef, 5*(ny#=*-2*=#)*(nz-2))
      buf_exch_n_in[z] = Array{FloatType}(undef, 5*nx*(nz-2))
      buf_exch_s_in[z] = Array{FloatType}(undef, 5*nx*(nz-2))

      buf_exch_w_out[z] = Array{FloatType}(undef, 5*(ny#=*-2*=#)*(nz-2))
      buf_exch_e_out[z] = Array{FloatType}(undef, 5*(ny#=*-2*=#)*(nz-2))
      buf_exch_n_out[z] = Array{FloatType}(undef, 5*nx*(nz-2))
      buf_exch_s_out[z] = Array{FloatType}(undef, 5*nx*(nz-2))

      return nothing
      
end

function alloc_field_space_zones(proc_num_zones)

      global rsdnm = Array{Array{FloatType}}(undef, proc_num_zones)
      global errnm = Array{Array{FloatType}}(undef, proc_num_zones)

      global u = Array{OffsetArray{FloatType, 4, Array{FloatType, 4}}}(undef, proc_num_zones)
      global rsd = Array{OffsetArray{FloatType, 4, Array{FloatType, 4}}}(undef, proc_num_zones)
      global frct = Array{OffsetArray{FloatType, 4, Array{FloatType, 4}}}(undef, proc_num_zones)
      global flux = Array{OffsetArray{FloatType, 4, Array{FloatType, 4}}}(undef, proc_num_zones)

      global a = Array{Array{FloatType, 3}}(undef, proc_num_zones)
      global b = Array{Array{FloatType, 3}}(undef, proc_num_zones) 
      global c = Array{Array{FloatType, 3}}(undef, proc_num_zones) 
      global d = Array{Array{FloatType, 3}}(undef, proc_num_zones) 

      global phi1 = Array{OffsetArray{FloatType, 2, Array{FloatType, 2}}}(undef, proc_num_zones)
      global phi2 = Array{OffsetArray{FloatType, 2, Array{FloatType, 2}}}(undef, proc_num_zones)

      global buf = Array{Array{FloatType, 2}}(undef, proc_num_zones) 
      global buf1 = Array{Array{FloatType, 2}}(undef, proc_num_zones) 

      global buf_exch_w_in = Array{Array{FloatType}}(undef, proc_num_zones) 
      global buf_exch_e_in = Array{Array{FloatType}}(undef, proc_num_zones) 
      global buf_exch_n_in = Array{Array{FloatType}}(undef, proc_num_zones) 
      global buf_exch_s_in = Array{Array{FloatType}}(undef, proc_num_zones) 

      global buf_exch_w_out = Array{Array{FloatType}}(undef, proc_num_zones) 
      global buf_exch_e_out = Array{Array{FloatType}}(undef, proc_num_zones) 
      global buf_exch_n_out = Array{Array{FloatType}}(undef, proc_num_zones) 
      global buf_exch_s_out = Array{Array{FloatType}}(undef, proc_num_zones) 

      return nothing
            
end