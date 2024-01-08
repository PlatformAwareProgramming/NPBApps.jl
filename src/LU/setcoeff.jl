#---------------------------------------------------------------------
#   set up coefficients
#---------------------------------------------------------------------


#---------------------------------------------------------------------
#   diffusion coefficients
#---------------------------------------------------------------------
const dx1 = 0.75e+00
const dx2 = dx1
const dx3 = dx1
const dx4 = dx1
const dx5 = dx1

const dy1 = 0.75e+00
const dy2 = dy1
const dy3 = dy1
const dy4 = dy1
const dy5 = dy1

const dz1 = 1.00e+00
const dz2 = dz1
const dz3 = dz1
const dz4 = dz1
const dz5 = dz1

const ii1 = 2
const ji1 = 2
const ki1 = 3

#---------------------------------------------------------------------
#   fourth difference dissipation
#---------------------------------------------------------------------
const dssp = ( max(dx1, dy1, dz1 ) ) / 4.0e+00

const ce = SA_F64[2.0e0 0.0e0 0.0e0 4.0e0 5.0e0 3.0e0 0.5e0 0.02e0 0.01e0 0.03e0 0.5e0 0.4e0 0.3e0;  # coefficients of the exact solution to the first pde
                  1.0e0 0.0e0 0.0e0 0.0e0 1.0e0 2.0e0 3.0e0 0.01e0 0.03e0 0.02e0 0.4e0 0.3e0 0.5e0;  # coefficients of the exact solution to the second pde
                  2.0e0 2.0e0 0.0e0 0.0e0 0.0e0 2.0e0 3.0e0 0.04e0 0.03e0 0.05e0 0.3e0 0.5e0 0.4e0;  # coefficients of the exact solution to the third pde
                  2.0e0 2.0e0 0.0e0 0.0e0 0.0e0 2.0e0 3.0e0 0.03e0 0.05e0 0.04e0 0.2e0 0.1e0 0.3e0;  # coefficients of the exact solution to the fourth pde
                  5.0e0 4.0e0 3.0e0 2.0e0 0.1e0 0.4e0 0.3e0 0.05e0 0.04e0 0.03e0 0.1e0 0.3e0 0.2e0]  # coefficients of the exact solution to the fifth pde


 function setcoeff()

#---------------------------------------------------------------------
#  local variables
#---------------------------------------------------------------------

      global dxi = 1.0e+00 / ( nx0 - 1 )
      global deta = 1.0e+00 / ( ny0 - 1 )
      global dzeta = 1.0e+00 / ( nz0 - 1 )

      global tx1 = 1.0e+00 / ( dxi * dxi )
      global tx2 = 1.0e+00 / ( 2.0e+00 * dxi )
      global tx3 = 1.0e+00 / dxi

      global ty1 = 1.0e+00 / ( deta * deta )
      global ty2 = 1.0e+00 / ( 2.0e+00 * deta )
      global ty3 = 1.0e+00 / deta

      global tz1 = 1.0e+00 / ( dzeta * dzeta )
      global tz2 = 1.0e+00 / ( 2.0e+00 * dzeta )
      global tz3 = 1.0e+00 / dzeta

      global ii2 = nx0 - 1
      global ji2 = ny0 - 2
      global ki2 = nz0 - 1

      return nothing
end


