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
const t_last = 10


#---------------------------------------------------------------------
# allocate space dynamically for data arrays
#---------------------------------------------------------------------

 function alloc_space()

      global rsdnm = Array{Float64}(undef,5)
      global errnm = Array{Float64}(undef,5)

      u0 = zeros(Float64, 5, isiz1 + 4, isiz2 + 4, isiz3)
      global u = OffsetArray(u0, 1:5, -1:isiz1+2, -1:isiz2+2, 1:isiz3)

      rsd0 = zeros(Float64, 5, isiz1 + 4, isiz2 + 4, isiz3)
      global rsd = OffsetArray(rsd0, 1:5, -1:isiz1+2, -1:isiz2+2, 1:isiz3)

      frct0 = zeros(Float64, 5, isiz1 + 4, isiz2 + 4, isiz3)
      global frct = OffsetArray(frct0, 1:5, -1:isiz1+2, -1:isiz2+2, 1:isiz3)

      flux0 = zeros(Float64, 5, isiz1 + 2, isiz2 + 2, isiz3)
      global flux = OffsetArray(flux0, 1:5, 0:isiz1+1,  0:isiz2+1, 1:isiz3)

      global a = zeros(Float64, 5, 5, isiz1)
      global b = zeros(Float64, 5, 5, isiz1)
      global c = zeros(Float64, 5, 5, isiz1)
      global d = zeros(Float64, 5, 5, isiz1)

      phi10 = zeros(Float64, isiz2 + 2, isiz3 + 2)
      global phi1 = OffsetArray(phi10, 0:isiz2+1, 0:isiz3+1)

      phi20 = zeros(Float64, isiz2 + 2, isiz3 + 2)
      global phi2 = OffsetArray(phi20, 0:isiz2+1, 0:isiz3+1)

      global buf = zeros(Float64, 5, 2*isiz2*isiz3)
      global buf1 = zeros(Float64, 5, 2*isiz2*isiz3)

      return nothing
      
end

