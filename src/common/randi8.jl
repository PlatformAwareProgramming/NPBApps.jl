
#---------------------------------------------------------------------
#
#   This routine returns a uniform pseudorandom double precision number in the
#   range (0, 1) by using the linear congruential generator
#
#   x_{k+1} = a x_k  (mod 2^46)
#
#   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
#   before repeating.  The argument A is the same as 'a' in the above formula,
#   and X is the same as x_0.  A and X must be odd double precision integers
#   in the range (1, 2^46).  The returned value RANDLC is normalized to be
#   between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain
#   the new seed x_1, so that subsequent calls to RANDLC using the same
#   arguments will generate a continuous sequence.

function randlc(x, a)

      d2m46 = 0.5e0^46

      i246m1 = truncComplex(Int, Z"00003FFFFFFFFFFF", 8)

      Lx = X
      La = A

      Lx   = iand(Lx*La, i246m1)
      randlc = d2m46*float(Lx)
      x    = float(Lx)
      return nothing
end


#---------------------------------------------------------------------
#---------------------------------------------------------------------

function VRANLC(N, X, A, Y)

# This doesn't work, because the compiler does the calculation in 32
# bits and overflows. No standard way (without f90 stuff) to specify
# that the rhs should be done in 64 bit arithmetic. 
#      parameter(i246m1=2**46-1)

      d2m46 = 0.5e0^46

      i246m1 = truncComplex(Int, Z"00003FFFFFFFFFFF", 8)

      Lx = X
      La = A
      for i = 1:N
         Lx   = iand(Lx*La, i246m1)
         y[i] = d2m46*float(Lx)
      end
      x    = float(Lx)

      return nothing
end

