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
#
#   This routine should produce the same results on any computer with at least
#   48 mantissa bits in double precision floating point data.  On 64 bit
#   systems, double precision should be disabled.
#
#   David H. Bailey     October 26, 1990
#
#---------------------------------------------------------------------

function randlc(x, a)

      r23 = 0.5e0 ^ 23; r46 = r23 ^ 2; t23 = 2.0e0 ^ 23;
      t46 = t23 ^ 2

#---------------------------------------------------------------------
#   Break A into two parts such that A = 2^23 * A1 + A2.
#---------------------------------------------------------------------
      t1 = r23 * a
      a1 = trunc(Int, t1)
      a2 = a - t23 * a1

#---------------------------------------------------------------------
#   Break X into two parts such that X = 2^23 * X1 + X2, compute
#   Z = A1 * X2 + A2 * X1  (mod 2^23), and then
#   X = 2^23 * Z + A2 * X2  (mod 2^46).
#---------------------------------------------------------------------
      t1 = r23 * x
      x1 = trunc(Int, t1)
      x2 = x - t23 * x1


      t1 = a1 * x2 + a2 * x1
      t2 = trunc(Int, r23 * t1)
      z = t1 - t23 * t2
      t3 = t23 * z + a2 * x2
      t4 = trunc(Int, r46 * t3)
      x = t3 - t46 * t4
      randlc = r46 * x
      return nothing
end


#---------------------------------------------------------------------
#   This routine generates N uniform pseudorandom double precision numbers in
#   the range (0, 1) by using the linear congruential generator
#   
#   x_{k+1} = a x_k  (mod 2^46)
#   
#   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
#   before repeating.  The argument A is the same as 'a' in the above formula,
#   and X is the same as x_0.  A and X must be odd double precision integers
#   in the range (1, 2^46).  The N results are placed in Y and are normalized
#   to be between 0 and 1.  X is updated to contain the new seed, so that
#   subsequent calls to RANDLC using the same arguments will generate a
#   continuous sequence.
#   
#   This routine generates the output sequence in batches of length NV, for
#   convenience on vector computers.  This routine should produce the same
#   results on any computer with at least 48 mantissa bits in double precision
#   floating point data.  On Cray systems, double precision should be disabled.
#   
#   David H. Bailey    August 30, 1990
#---------------------------------------------------------------------

function vranlc(n, x, a, y)

      r23 = 2.0e0 ^ Complex(-23); r46 = r23 * r23; t23 = 2.0e0 ^ 23;
      t46 = t23 * t23; nv = 64

#---------------------------------------------------------------------
#     Compute the first NV elements of the sequence using RANDLC.
#---------------------------------------------------------------------
      t1 = x
      n1 = min(n, nv)

      for i = 1:n1
         xv[i] = t46 * randlc(t1, a)
      end

#---------------------------------------------------------------------
#     It is not necessary to compute AN, A1 or A2 unless N is greater than NV.
#---------------------------------------------------------------------
      if n > nv

#---------------------------------------------------------------------
#     Compute AN = AA ^ NV (mod 2^46) using successive calls to RANDLC.
#---------------------------------------------------------------------
         t1 = a
         t2 = r46 * a

         for i = 1:nv - 1
            t2 = randlc(t1, a)
         end

         an = t46 * t2

#---------------------------------------------------------------------
#     Break AN into two parts such that AN = 2^23 * A1 + A2.
#---------------------------------------------------------------------
         t1 = r23 * an
         a1 = aint(t1)
         a2 = an - t23 * a1
      end

#---------------------------------------------------------------------
#     Compute N pseudorandom results in batches of size NV.
#---------------------------------------------------------------------
      for j = 0:nv:n - 1
         n1 = min(nv, n - j)

#---------------------------------------------------------------------
#     Compute up to NV results based on the current seed vector XV.
#---------------------------------------------------------------------
         for i = 1:n1
            y[i+j] = r46 * xv[i]
         end

#---------------------------------------------------------------------
#     If this is the last pass through the 140 loop, it is not necessary to
#     update the XV vector.
#---------------------------------------------------------------------
         if (j + n1 == n) @goto L150 end

#---------------------------------------------------------------------
#     Update the XV vector by multiplying each element by AN (mod 2^46).
#---------------------------------------------------------------------
         for i = 1:nv
            t1 = r23 * xv[i]
            x1 = aint(t1)
            x2 = xv[i] - t23 * x1
            t1 = a1 * x2 + a2 * x1
            t2 = aint(r23 * t1)
            yy = t1 - t23 * t2
            t3 = t23 * yy + a2 * x2
            t4 = aint(r46 * t3)
            xv[i] = t3 - t46 * t4
         end

      end

#---------------------------------------------------------------------
#     Save the last seed in X so that subsequent calls to VRANLC will generate
#     a continuous sequence.
#---------------------------------------------------------------------
      @label L150
      x = xv[n1]

      return nothing
end

#----- end of program ------------------------------------------------

