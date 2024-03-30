#---------------------------------------------------------------------
#
#   compute the regular-sparse, block upper triangular solution:
#
#                     v <-- ( U-inv ) * v
#
#---------------------------------------------------------------------

 function buts( j, k, omega, v, tv, d, udx, udy, udz, ist, iend )

#---------------------------------------------------------------------
#---------------------------------------------------------------------

         for i = iend:-1:ist
            for m = 1:5
                  tv[m,i] = omega * (  udz[m,1,i] * v[1,i,j,k+1]+
                                       udz[m,2,i] * v[2,i,j,k+1]+
                                       udz[m,3,i] * v[3,i,j,k+1]+
                                       udz[m,4,i] * v[4,i,j,k+1]+
                                       udz[m,5,i] * v[5,i,j,k+1] )
            end
         end


         for i = iend:-1:ist

            for m = 1:5
                  tv[m,i] += omega * ( udy[m,1,i] * v[1,i,j+1,k]+
                                       udx[m,1,i] * v[1,i+1,j,k]+
                                       udy[m,2,i] * v[2,i,j+1,k]+
                                       udx[m,2,i] * v[2,i+1,j,k]+
                                       udy[m,3,i] * v[3,i,j+1,k]+
                                       udx[m,3,i] * v[3,i+1,j,k]+
                                       udy[m,4,i] * v[4,i,j+1,k]+
                                       udx[m,4,i] * v[4,i+1,j,k]+
                                       udy[m,5,i] * v[5,i,j+1,k]+
                                       udx[m,5,i] * v[5,i+1,j,k] )
            end

#---------------------------------------------------------------------
#   diagonal block inversion
#---------------------------------------------------------------------
            for m = 1:5
               tmat[m,1] = d[m,1,i]
               tmat[m,2] = d[m,2,i]
               tmat[m,3] = d[m,3,i]
               tmat[m,4] = d[m,4,i]
               tmat[m,5] = d[m,5,i]
            end

            tmp1 = 1.0e+00 / tmat[1,1]
            tmp = tmp1 * tmat[2,1]
            tmat[2,2] -= tmp * tmat[1,2]
            tmat[2,3] -= tmp * tmat[1,3]
            tmat[2,4] -= tmp * tmat[1,4]
            tmat[2,5] -= tmp * tmat[1,5]
            tv[2,i] -= tv[1,i] * tmp

            tmp = tmp1 * tmat[3,1]
            tmat[3,2] -= tmp * tmat[1,2]
            tmat[3,3] -= tmp * tmat[1,3]
            tmat[3,4] -= tmp * tmat[1,4]
            tmat[3,5] -= tmp * tmat[1,5]
            tv[3,i] -= tv[1,i] * tmp

            tmp = tmp1 * tmat[4,1]
            tmat[4,2] -= tmp * tmat[1,2]
            tmat[4,3] -= tmp * tmat[1,3]
            tmat[4,4] -= tmp * tmat[1,4]
            tmat[4,5] -= tmp * tmat[1,5]
            tv[4,i] -= tv[1,i] * tmp

            tmp = tmp1 * tmat[5,1]
            tmat[5,2] -= tmp * tmat[1,2]
            tmat[5,3] -= tmp * tmat[1,3]
            tmat[5,4] -= tmp * tmat[1,4]
            tmat[5,5] -= tmp * tmat[1,5]
            tv[5,i] -= tv[1,i] * tmp


            tmp1 = 1.0e+00 / tmat[2,2]
            tmp = tmp1 * tmat[3,2]
            tmat[3,3] -= tmp * tmat[2,3]
            tmat[3,4] -= tmp * tmat[2,4]
            tmat[3,5] -= tmp * tmat[2,5]
            tv[3,i] = tv[3,i] - tv[2,i] * tmp

            tmp = tmp1 * tmat[4,2]
            tmat[4,3] -= tmp * tmat[2,3]
            tmat[4,4] -= tmp * tmat[2,4]
            tmat[4,5] -= tmp * tmat[2,5]
            tv[4,i] -= tv[2,i] * tmp

            tmp = tmp1 * tmat[5,2]
            tmat[5,3] -= tmp * tmat[2,3]
            tmat[5,4] -= tmp * tmat[2,4]
            tmat[5,5] -=- tmp * tmat[2,5]
            tv[5,i] -= tv[2,i] * tmp


            tmp1 = 1.0e+00 / tmat[3,3]
            tmp = tmp1 * tmat[4,3]
            tmat[4,4] -= tmp * tmat[3,4]
            tmat[4,5] -= tmp * tmat[3,5]
            tv[4,i] -= tv[3,i] * tmp

            tmp = tmp1 * tmat[5,3]
            tmat[5,4] -=tmp * tmat[3,4]
            tmat[5,5] -= tmp * tmat[3,5]
            tv[5,i] -= tv[3,i] * tmp


            tmp1 = 1.0e+00 / tmat[4,4]
            tmp = tmp1 * tmat[5,4]
            tmat[5,5] -= tmp * tmat[4,5]
            tv[5,i] -= tv[4,i] * tmp

#---------------------------------------------------------------------
#   back substitution
#---------------------------------------------------------------------
            tv[5,i] /= tmat[5,5]

            tv[4,i] -= tmat[4,5] * tv[5,i]
            tv[4,i] /= tmat[4,4]

            tv[3,i] -= tmat[3,4] * tv[4,i] - tmat[3,5] * tv[5,i]
            tv[3,i] /= tmat[3,3]

            tv[2,i] -= tmat[2,3] * tv[3,i] - tmat[2,4] * tv[4,i] - tmat[2,5] * tv[5,i]
            tv[2,i] /= tmat[2,2]

            tv[1,i] -= tmat[1,2] * tv[2,i] - tmat[1,3] * tv[3,i] - tmat[1,4] * tv[4,i] - tmat[1,5] * tv[5,i]
            tv[1,i] /= tmat[1,1]

            v[1,i,j,k] -= tv[1,i]
            v[2,i,j,k] -= tv[2,i]
            v[3,i,j,k] -= tv[3,i]
            v[4,i,j,k] -= tv[4,i]
            v[5,i,j,k] -= tv[5,i]


         end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

      return nothing
end
