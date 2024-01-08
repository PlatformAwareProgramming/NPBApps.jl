#---------------------------------------------------------------------
#     subtracts bvec=bvec - ablock*avec
#---------------------------------------------------------------------
function matvec_sub(ablock, i, avec, a1, a2, a3, a4, bvec, b1, b2, b3, b4)


#---------------------------------------------------------------------
#            rhs[i,ic,jc,kc,ccell] = rhs[i,ic,jc,kc,ccell] 
#     $           - lhs(i,1,ablock,ia,ja,ka,acell)*
#---------------------------------------------------------------------
       bvec[1,b1,b2,b3,b4] -= ablock[1,1,i]*avec[1,a1,a2,a3,a4] +
                              ablock[1,2,i]*avec[2,a1,a2,a3,a4] +
                              ablock[1,3,i]*avec[3,a1,a2,a3,a4] +
                              ablock[1,4,i]*avec[4,a1,a2,a3,a4] +
                              ablock[1,5,i]*avec[5,a1,a2,a3,a4]
       bvec[2,b1,b2,b3,b4] -= ablock[2,1,i]*avec[1,a1,a2,a3,a4] +
                              ablock[2,2,i]*avec[2,a1,a2,a3,a4] +
                              ablock[2,3,i]*avec[3,a1,a2,a3,a4] +
                              ablock[2,4,i]*avec[4,a1,a2,a3,a4] +
                              ablock[2,5,i]*avec[5,a1,a2,a3,a4]
       bvec[3,b1,b2,b3,b4] -= ablock[3,1,i]*avec[1,a1,a2,a3,a4] +
                              ablock[3,2,i]*avec[2,a1,a2,a3,a4] +
                              ablock[3,3,i]*avec[3,a1,a2,a3,a4] +
                              ablock[3,4,i]*avec[4,a1,a2,a3,a4] +
                              ablock[3,5,i]*avec[5,a1,a2,a3,a4]
       bvec[4,b1,b2,b3,b4] -= ablock[4,1,i]*avec[1,a1,a2,a3,a4] +
                              ablock[4,2,i]*avec[2,a1,a2,a3,a4] +
                              ablock[4,3,i]*avec[3,a1,a2,a3,a4] +
                              ablock[4,4,i]*avec[4,a1,a2,a3,a4] +
                              ablock[4,5,i]*avec[5,a1,a2,a3,a4]
       bvec[5,b1,b2,b3,b4] -= ablock[5,1,i]*avec[1,a1,a2,a3,a4] +
                              ablock[5,2,i]*avec[2,a1,a2,a3,a4] +
                              ablock[5,3,i]*avec[3,a1,a2,a3,a4] +
                              ablock[5,4,i]*avec[4,a1,a2,a3,a4] +
                              ablock[5,5,i]*avec[5,a1,a2,a3,a4]
      return nothing
end

#---------------------------------------------------------------------
#     subtracts a[i,j,k] X  b[i,j,k] from c[i,j,k]
#---------------------------------------------------------------------

function matmul_sub(ablock, ia, bblock, b1, b2, b3, b4, cblock, ic)

      timer_start(99)

     #= for j = 1:5, i = 1:5
                  sum = 0
                  for k = 1:5
                      sum += ablock[i,k,ia]*bblock[k,j,b1,b2,b3,b4]
                  end
                  cblock[i,j,ic] -= sum
            #end
      end=#

      
      cblock[1,1,ic] -= ablock[1,1,ia]*bblock[1,1,b1,b2,b3,b4] +
                        ablock[1,2,ia]*bblock[2,1,b1,b2,b3,b4] +
                        ablock[1,3,ia]*bblock[3,1,b1,b2,b3,b4] +
                        ablock[1,4,ia]*bblock[4,1,b1,b2,b3,b4] +
                        ablock[1,5,ia]*bblock[5,1,b1,b2,b3,b4]
      cblock[2,1,ic] -= ablock[2,1,ia]*bblock[1,1,b1,b2,b3,b4] +
                        ablock[2,2,ia]*bblock[2,1,b1,b2,b3,b4] +
                        ablock[2,3,ia]*bblock[3,1,b1,b2,b3,b4] +
                        ablock[2,4,ia]*bblock[4,1,b1,b2,b3,b4] +
                        ablock[2,5,ia]*bblock[5,1,b1,b2,b3,b4]
      cblock[3,1,ic] -= ablock[3,1,ia]*bblock[1,1,b1,b2,b3,b4] +
                        ablock[3,2,ia]*bblock[2,1,b1,b2,b3,b4] +
                        ablock[3,3,ia]*bblock[3,1,b1,b2,b3,b4] +
                        ablock[3,4,ia]*bblock[4,1,b1,b2,b3,b4] +
                        ablock[3,5,ia]*bblock[5,1,b1,b2,b3,b4]
      cblock[4,1,ic] -= ablock[4,1,ia]*bblock[1,1,b1,b2,b3,b4] +
                        ablock[4,2,ia]*bblock[2,1,b1,b2,b3,b4] +
                        ablock[4,3,ia]*bblock[3,1,b1,b2,b3,b4] +
                        ablock[4,4,ia]*bblock[4,1,b1,b2,b3,b4] +
                        ablock[4,5,ia]*bblock[5,1,b1,b2,b3,b4]
      cblock[5,1,ic] -= ablock[5,1,ia]*bblock[1,1,b1,b2,b3,b4] +
                        ablock[5,2,ia]*bblock[2,1,b1,b2,b3,b4] +
                        ablock[5,3,ia]*bblock[3,1,b1,b2,b3,b4] +
                        ablock[5,4,ia]*bblock[4,1,b1,b2,b3,b4] +
                        ablock[5,5,ia]*bblock[5,1,b1,b2,b3,b4]
      cblock[1,2,ic] -= ablock[1,1,ia]*bblock[1,2,b1,b2,b3,b4] +
                        ablock[1,2,ia]*bblock[2,2,b1,b2,b3,b4] +
                        ablock[1,3,ia]*bblock[3,2,b1,b2,b3,b4] +
                        ablock[1,4,ia]*bblock[4,2,b1,b2,b3,b4] +
                        ablock[1,5,ia]*bblock[5,2,b1,b2,b3,b4]
      cblock[2,2,ic] -= ablock[2,1,ia]*bblock[1,2,b1,b2,b3,b4] +
                        ablock[2,2,ia]*bblock[2,2,b1,b2,b3,b4] +
                        ablock[2,3,ia]*bblock[3,2,b1,b2,b3,b4] +
                        ablock[2,4,ia]*bblock[4,2,b1,b2,b3,b4] +
                        ablock[2,5,ia]*bblock[5,2,b1,b2,b3,b4]
      cblock[3,2,ic] -= ablock[3,1,ia]*bblock[1,2,b1,b2,b3,b4] +
                        ablock[3,2,ia]*bblock[2,2,b1,b2,b3,b4] +
                        ablock[3,3,ia]*bblock[3,2,b1,b2,b3,b4] +
                        ablock[3,4,ia]*bblock[4,2,b1,b2,b3,b4] +
                        ablock[3,5,ia]*bblock[5,2,b1,b2,b3,b4]
      cblock[4,2,ic] -= ablock[4,1,ia]*bblock[1,2,b1,b2,b3,b4] +
                        ablock[4,2,ia]*bblock[2,2,b1,b2,b3,b4] +
                        ablock[4,3,ia]*bblock[3,2,b1,b2,b3,b4] +
                        ablock[4,4,ia]*bblock[4,2,b1,b2,b3,b4] +
                        ablock[4,5,ia]*bblock[5,2,b1,b2,b3,b4]
      cblock[5,2,ic] -= ablock[5,1,ia]*bblock[1,2,b1,b2,b3,b4] +
                        ablock[5,2,ia]*bblock[2,2,b1,b2,b3,b4] +
                        ablock[5,3,ia]*bblock[3,2,b1,b2,b3,b4] +
                        ablock[5,4,ia]*bblock[4,2,b1,b2,b3,b4] +
                        ablock[5,5,ia]*bblock[5,2,b1,b2,b3,b4]
      cblock[1,3,ic] -= ablock[1,1,ia]*bblock[1,3,b1,b2,b3,b4] +
                        ablock[1,2,ia]*bblock[2,3,b1,b2,b3,b4] +
                        ablock[1,3,ia]*bblock[3,3,b1,b2,b3,b4] +
                        ablock[1,4,ia]*bblock[4,3,b1,b2,b3,b4] +
                        ablock[1,5,ia]*bblock[5,3,b1,b2,b3,b4]
      cblock[2,3,ic] -= ablock[2,1,ia]*bblock[1,3,b1,b2,b3,b4] +
                        ablock[2,2,ia]*bblock[2,3,b1,b2,b3,b4] +
                        ablock[2,3,ia]*bblock[3,3,b1,b2,b3,b4] +
                        ablock[2,4,ia]*bblock[4,3,b1,b2,b3,b4] +
                        ablock[2,5,ia]*bblock[5,3,b1,b2,b3,b4]
      cblock[3,3,ic] -= ablock[3,1,ia]*bblock[1,3,b1,b2,b3,b4] +
                        ablock[3,2,ia]*bblock[2,3,b1,b2,b3,b4] +
                        ablock[3,3,ia]*bblock[3,3,b1,b2,b3,b4] +
                        ablock[3,4,ia]*bblock[4,3,b1,b2,b3,b4] +
                        ablock[3,5,ia]*bblock[5,3,b1,b2,b3,b4]
      cblock[4,3,ic] -= ablock[4,1,ia]*bblock[1,3,b1,b2,b3,b4] +
                        ablock[4,2,ia]*bblock[2,3,b1,b2,b3,b4] +
                        ablock[4,3,ia]*bblock[3,3,b1,b2,b3,b4] +
                        ablock[4,4,ia]*bblock[4,3,b1,b2,b3,b4] +
                        ablock[4,5,ia]*bblock[5,3,b1,b2,b3,b4]
      cblock[5,3,ic] -= ablock[5,1,ia]*bblock[1,3,b1,b2,b3,b4] +
                        ablock[5,2,ia]*bblock[2,3,b1,b2,b3,b4] +
                        ablock[5,3,ia]*bblock[3,3,b1,b2,b3,b4] +
                        ablock[5,4,ia]*bblock[4,3,b1,b2,b3,b4] +
                        ablock[5,5,ia]*bblock[5,3,b1,b2,b3,b4]
      cblock[1,4,ic] -= ablock[1,1,ia]*bblock[1,4,b1,b2,b3,b4] +
                        ablock[1,2,ia]*bblock[2,4,b1,b2,b3,b4] +
                        ablock[1,3,ia]*bblock[3,4,b1,b2,b3,b4] +
                        ablock[1,4,ia]*bblock[4,4,b1,b2,b3,b4] +
                        ablock[1,5,ia]*bblock[5,4,b1,b2,b3,b4]
      cblock[2,4,ic] -= ablock[2,1,ia]*bblock[1,4,b1,b2,b3,b4] +
                        ablock[2,2,ia]*bblock[2,4,b1,b2,b3,b4] +
                        ablock[2,3,ia]*bblock[3,4,b1,b2,b3,b4] +
                        ablock[2,4,ia]*bblock[4,4,b1,b2,b3,b4] +
                        ablock[2,5,ia]*bblock[5,4,b1,b2,b3,b4]
      cblock[3,4,ic] -= ablock[3,1,ia]*bblock[1,4,b1,b2,b3,b4] +
                        ablock[3,2,ia]*bblock[2,4,b1,b2,b3,b4] +
                        ablock[3,3,ia]*bblock[3,4,b1,b2,b3,b4] +
                        ablock[3,4,ia]*bblock[4,4,b1,b2,b3,b4] +
                        ablock[3,5,ia]*bblock[5,4,b1,b2,b3,b4]
      cblock[4,4,ic] -= ablock[4,1,ia]*bblock[1,4,b1,b2,b3,b4] +
                        ablock[4,2,ia]*bblock[2,4,b1,b2,b3,b4] +
                        ablock[4,3,ia]*bblock[3,4,b1,b2,b3,b4] +
                        ablock[4,4,ia]*bblock[4,4,b1,b2,b3,b4] +
                        ablock[4,5,ia]*bblock[5,4,b1,b2,b3,b4]
      cblock[5,4,ic] -= ablock[5,1,ia]*bblock[1,4,b1,b2,b3,b4] +
                        ablock[5,2,ia]*bblock[2,4,b1,b2,b3,b4] +
                        ablock[5,3,ia]*bblock[3,4,b1,b2,b3,b4] +
                        ablock[5,4,ia]*bblock[4,4,b1,b2,b3,b4] +
                        ablock[5,5,ia]*bblock[5,4,b1,b2,b3,b4]
      cblock[1,5,ic] -= ablock[1,1,ia]*bblock[1,5,b1,b2,b3,b4] +
                        ablock[1,2,ia]*bblock[2,5,b1,b2,b3,b4] +
                        ablock[1,3,ia]*bblock[3,5,b1,b2,b3,b4] +
                        ablock[1,4,ia]*bblock[4,5,b1,b2,b3,b4] +
                        ablock[1,5,ia]*bblock[5,5,b1,b2,b3,b4]
      cblock[2,5,ic] -= ablock[2,1,ia]*bblock[1,5,b1,b2,b3,b4] +
                        ablock[2,2,ia]*bblock[2,5,b1,b2,b3,b4] +
                        ablock[2,3,ia]*bblock[3,5,b1,b2,b3,b4] +
                        ablock[2,4,ia]*bblock[4,5,b1,b2,b3,b4] +
                        ablock[2,5,ia]*bblock[5,5,b1,b2,b3,b4]
      cblock[3,5,ic] -= ablock[3,1,ia]*bblock[1,5,b1,b2,b3,b4] +
                        ablock[3,2,ia]*bblock[2,5,b1,b2,b3,b4] +
                        ablock[3,3,ia]*bblock[3,5,b1,b2,b3,b4] +
                        ablock[3,4,ia]*bblock[4,5,b1,b2,b3,b4] +
                        ablock[3,5,ia]*bblock[5,5,b1,b2,b3,b4]
      cblock[4,5,ic] -= ablock[4,1,ia]*bblock[1,5,b1,b2,b3,b4] +
                        ablock[4,2,ia]*bblock[2,5,b1,b2,b3,b4] +
                        ablock[4,3,ia]*bblock[3,5,b1,b2,b3,b4] +
                        ablock[4,4,ia]*bblock[4,5,b1,b2,b3,b4] +
                        ablock[4,5,ia]*bblock[5,5,b1,b2,b3,b4]
      cblock[5,5,ic] -= ablock[5,1,ia]*bblock[1,5,b1,b2,b3,b4] +
                        ablock[5,2,ia]*bblock[2,5,b1,b2,b3,b4] +
                        ablock[5,3,ia]*bblock[3,5,b1,b2,b3,b4] +
                        ablock[5,4,ia]*bblock[4,5,b1,b2,b3,b4] +
                        ablock[5,5,ia]*bblock[5,5,b1,b2,b3,b4]


      timer_stop(99)

      return nothing
end 

#---------------------------------------------------------------------
#     subtracts a[i,j,k] X  b[i,j,k] from c[i,j,k]
#---------------------------------------------------------------------

function matmul_sub2(ablock, bblock, cblock)


      timer_start(99)

#=     for j = 1:5, i = 1:5
                  sum = 0
                  for k = 1:5
                      sum += ablock[i,k]*bblock[k,j]
                  end
                  cblock[i,j] -= sum
            #end
      end
=#
         cblock[1,1] -= ablock[1,1]*bblock[1,1] +
                        ablock[1,2]*bblock[2,1] +
                        ablock[1,3]*bblock[3,1] +
                        ablock[1,4]*bblock[4,1] +
                        ablock[1,5]*bblock[5,1]
         cblock[2,1] -= ablock[2,1]*bblock[1,1] +
                        ablock[2,2]*bblock[2,1] +
                        ablock[2,3]*bblock[3,1] +
                        ablock[2,4]*bblock[4,1] +
                        ablock[2,5]*bblock[5,1]
         cblock[3,1] -= ablock[3,1]*bblock[1,1] +
                        ablock[3,2]*bblock[2,1] +
                        ablock[3,3]*bblock[3,1] +
                        ablock[3,4]*bblock[4,1] +
                        ablock[3,5]*bblock[5,1]
         cblock[4,1] -= ablock[4,1]*bblock[1,1] +
                        ablock[4,2]*bblock[2,1] +
                        ablock[4,3]*bblock[3,1] +
                        ablock[4,4]*bblock[4,1] +
                        ablock[4,5]*bblock[5,1]
         cblock[5,1] -= ablock[5,1]*bblock[1,1] +
                        ablock[5,2]*bblock[2,1] +
                        ablock[5,3]*bblock[3,1] +
                        ablock[5,4]*bblock[4,1] +
                        ablock[5,5]*bblock[5,1]
         cblock[1,2] -= ablock[1,1]*bblock[1,2] +
                        ablock[1,2]*bblock[2,2] +
                        ablock[1,3]*bblock[3,2] +
                        ablock[1,4]*bblock[4,2] +
                        ablock[1,5]*bblock[5,2]
         cblock[2,2] -= ablock[2,1]*bblock[1,2] +
                        ablock[2,2]*bblock[2,2] +
                        ablock[2,3]*bblock[3,2] +
                        ablock[2,4]*bblock[4,2] +
                        ablock[2,5]*bblock[5,2]
         cblock[3,2] -= ablock[3,1]*bblock[1,2] +
                        ablock[3,2]*bblock[2,2] +
                        ablock[3,3]*bblock[3,2] +
                        ablock[3,4]*bblock[4,2] +
                        ablock[3,5]*bblock[5,2]
         cblock[4,2] -= ablock[4,1]*bblock[1,2] +
                        ablock[4,2]*bblock[2,2] +
                        ablock[4,3]*bblock[3,2] +
                        ablock[4,4]*bblock[4,2] +
                        ablock[4,5]*bblock[5,2]
         cblock[5,2] -= ablock[5,1]*bblock[1,2] +
                        ablock[5,2]*bblock[2,2] +
                        ablock[5,3]*bblock[3,2] +
                        ablock[5,4]*bblock[4,2] +
                        ablock[5,5]*bblock[5,2]
         cblock[1,3] -= ablock[1,1]*bblock[1,3] +
                        ablock[1,2]*bblock[2,3] +
                        ablock[1,3]*bblock[3,3] +
                        ablock[1,4]*bblock[4,3] +
                        ablock[1,5]*bblock[5,3]
         cblock[2,3] -= ablock[2,1]*bblock[1,3] +
                        ablock[2,2]*bblock[2,3] +
                        ablock[2,3]*bblock[3,3] +
                        ablock[2,4]*bblock[4,3] +
                        ablock[2,5]*bblock[5,3]
         cblock[3,3] -= ablock[3,1]*bblock[1,3] +
                        ablock[3,2]*bblock[2,3] +
                        ablock[3,3]*bblock[3,3] +
                        ablock[3,4]*bblock[4,3] +
                        ablock[3,5]*bblock[5,3]
         cblock[4,3] -= ablock[4,1]*bblock[1,3] +
                        ablock[4,2]*bblock[2,3] +
                        ablock[4,3]*bblock[3,3] +
                        ablock[4,4]*bblock[4,3] +
                        ablock[4,5]*bblock[5,3]
         cblock[5,3] -= ablock[5,1]*bblock[1,3] +
                        ablock[5,2]*bblock[2,3] +
                        ablock[5,3]*bblock[3,3] +
                        ablock[5,4]*bblock[4,3] +
                        ablock[5,5]*bblock[5,3]
         cblock[1,4] -= ablock[1,1]*bblock[1,4] +
                        ablock[1,2]*bblock[2,4] +
                        ablock[1,3]*bblock[3,4] +
                        ablock[1,4]*bblock[4,4] +
                        ablock[1,5]*bblock[5,4]
         cblock[2,4] -= ablock[2,1]*bblock[1,4] +
                        ablock[2,2]*bblock[2,4] +
                        ablock[2,3]*bblock[3,4] +
                        ablock[2,4]*bblock[4,4] +
                        ablock[2,5]*bblock[5,4]
         cblock[3,4] -= ablock[3,1]*bblock[1,4] +
                        ablock[3,2]*bblock[2,4] +
                        ablock[3,3]*bblock[3,4] +
                        ablock[3,4]*bblock[4,4] +
                        ablock[3,5]*bblock[5,4]
         cblock[4,4] -= ablock[4,1]*bblock[1,4] +
                        ablock[4,2]*bblock[2,4] +
                        ablock[4,3]*bblock[3,4] +
                        ablock[4,4]*bblock[4,4] +
                        ablock[4,5]*bblock[5,4]
         cblock[5,4] -= ablock[5,1]*bblock[1,4] +
                        ablock[5,2]*bblock[2,4] +
                        ablock[5,3]*bblock[3,4] +
                        ablock[5,4]*bblock[4,4] +
                        ablock[5,5]*bblock[5,4]
         cblock[1,5] -= ablock[1,1]*bblock[1,5] +
                        ablock[1,2]*bblock[2,5] +
                        ablock[1,3]*bblock[3,5] +
                        ablock[1,4]*bblock[4,5] +
                        ablock[1,5]*bblock[5,5]
         cblock[2,5] -= ablock[2,1]*bblock[1,5] +
                        ablock[2,2]*bblock[2,5] +
                        ablock[2,3]*bblock[3,5] +
                        ablock[2,4]*bblock[4,5] +
                        ablock[2,5]*bblock[5,5]
         cblock[3,5] -= ablock[3,1]*bblock[1,5] +
                        ablock[3,2]*bblock[2,5] +
                        ablock[3,3]*bblock[3,5] +
                        ablock[3,4]*bblock[4,5] +
                        ablock[3,5]*bblock[5,5]
         cblock[4,5] -= ablock[4,1]*bblock[1,5] +
                        ablock[4,2]*bblock[2,5] +
                        ablock[4,3]*bblock[3,5] +
                        ablock[4,4]*bblock[4,5] +
                        ablock[4,5]*bblock[5,5]
         cblock[5,5] -= ablock[5,1]*bblock[1,5] +
                        ablock[5,2]*bblock[2,5] +
                        ablock[5,3]*bblock[3,5] +
                        ablock[5,4]*bblock[4,5] +
                        ablock[5,5]*bblock[5,5]

      timer_stop(99)

      return nothing
end



#---------------------------------------------------------------------
#---------------------------------------------------------------------

function binvcrhs(lhs, i, c, ic1, ic2, ic3, ic4, r, ir1, ir2, ir3, ir4)

      pivot = 1.00e0/lhs[1,1,i]
      lhs[1,2,i] *= pivot
      lhs[1,3,i] *= pivot
      lhs[1,4,i] *= pivot
      lhs[1,5,i] *= pivot
      c[1,1,ic1,ic2,ic3,ic4] *= pivot
      c[1,2,ic1,ic2,ic3,ic4] *= pivot
      c[1,3,ic1,ic2,ic3,ic4] *= pivot
      c[1,4,ic1,ic2,ic3,ic4] *= pivot
      c[1,5,ic1,ic2,ic3,ic4] *= pivot
      r[1,ir1,ir2,ir3,ir4]   *= pivot

      coeff = lhs[2,1,i]
      lhs[2,2,i] -= coeff*lhs[1,2,i]
      lhs[2,3,i] -= coeff*lhs[1,3,i]
      lhs[2,4,i] -= coeff*lhs[1,4,i]
      lhs[2,5,i] -= coeff*lhs[1,5,i]
      c[2,1,ic1,ic2,ic3,ic4] -= coeff*c[1,1,ic1,ic2,ic3,ic4]
      c[2,2,ic1,ic2,ic3,ic4] -= coeff*c[1,2,ic1,ic2,ic3,ic4]
      c[2,3,ic1,ic2,ic3,ic4] -= coeff*c[1,3,ic1,ic2,ic3,ic4]
      c[2,4,ic1,ic2,ic3,ic4] -= coeff*c[1,4,ic1,ic2,ic3,ic4]
      c[2,5,ic1,ic2,ic3,ic4] -= coeff*c[1,5,ic1,ic2,ic3,ic4]
      r[2,ir1,ir2,ir3,ir4] -= coeff*r[1,ir1,ir2,ir3,ir4]

      coeff = lhs[3,1,i]
      lhs[3,2,i] -= coeff*lhs[1,2,i]
      lhs[3,3,i] -= coeff*lhs[1,3,i]
      lhs[3,4,i] -= coeff*lhs[1,4,i]
      lhs[3,5,i] -= coeff*lhs[1,5,i]
      c[3,1,ic1,ic2,ic3,ic4] -= coeff*c[1,1,ic1,ic2,ic3,ic4]
      c[3,2,ic1,ic2,ic3,ic4] -= coeff*c[1,2,ic1,ic2,ic3,ic4]
      c[3,3,ic1,ic2,ic3,ic4] -= coeff*c[1,3,ic1,ic2,ic3,ic4]
      c[3,4,ic1,ic2,ic3,ic4] -= coeff*c[1,4,ic1,ic2,ic3,ic4]
      c[3,5,ic1,ic2,ic3,ic4] -= coeff*c[1,5,ic1,ic2,ic3,ic4]
      r[3,ir1,ir2,ir3,ir4] -= coeff*r[1,ir1,ir2,ir3,ir4]

      coeff = lhs[4,1,i]
      lhs[4,2,i] -= coeff*lhs[1,2,i]
      lhs[4,3,i] -= coeff*lhs[1,3,i]
      lhs[4,4,i] -= coeff*lhs[1,4,i]
      lhs[4,5,i] -= coeff*lhs[1,5,i]
      c[4,1,ic1,ic2,ic3,ic4] -= coeff*c[1,1,ic1,ic2,ic3,ic4]
      c[4,2,ic1,ic2,ic3,ic4] -= coeff*c[1,2,ic1,ic2,ic3,ic4]
      c[4,3,ic1,ic2,ic3,ic4] -= coeff*c[1,3,ic1,ic2,ic3,ic4]
      c[4,4,ic1,ic2,ic3,ic4] -= coeff*c[1,4,ic1,ic2,ic3,ic4]
      c[4,5,ic1,ic2,ic3,ic4] -= coeff*c[1,5,ic1,ic2,ic3,ic4]
      r[4,ir1,ir2,ir3,ir4] -= coeff*r[1,ir1,ir2,ir3,ir4]

      coeff = lhs[5,1,i]
      lhs[5,2,i] -= coeff*lhs[1,2,i]
      lhs[5,3,i] -= coeff*lhs[1,3,i]
      lhs[5,4,i] -= coeff*lhs[1,4,i]
      lhs[5,5,i] -= coeff*lhs[1,5,i]
      c[5,1,ic1,ic2,ic3,ic4] -= coeff*c[1,1,ic1,ic2,ic3,ic4]
      c[5,2,ic1,ic2,ic3,ic4] -= coeff*c[1,2,ic1,ic2,ic3,ic4]
      c[5,3,ic1,ic2,ic3,ic4] -= coeff*c[1,3,ic1,ic2,ic3,ic4]
      c[5,4,ic1,ic2,ic3,ic4] -= coeff*c[1,4,ic1,ic2,ic3,ic4]
      c[5,5,ic1,ic2,ic3,ic4] -= coeff*c[1,5,ic1,ic2,ic3,ic4]
      r[5,ir1,ir2,ir3,ir4] -= coeff*r[1,ir1,ir2,ir3,ir4]


      pivot = 1.00e0/lhs[2,2,i]
      lhs[2,3,i] *= pivot
      lhs[2,4,i] *= pivot
      lhs[2,5,i] *= pivot
      c[2,1,ic1,ic2,ic3,ic4] *= pivot
      c[2,2,ic1,ic2,ic3,ic4] *= pivot
      c[2,3,ic1,ic2,ic3,ic4] *= pivot
      c[2,4,ic1,ic2,ic3,ic4] *= pivot
      c[2,5,ic1,ic2,ic3,ic4] *= pivot
      r[2,ir1,ir2,ir3,ir4] *= pivot

      coeff = lhs[1,2,i]
      lhs[1,3,i] -= coeff*lhs[2,3,i]
      lhs[1,4,i] -= coeff*lhs[2,4,i]
      lhs[1,5,i] -= coeff*lhs[2,5,i]
      c[1,1,ic1,ic2,ic3,ic4] -= coeff*c[2,1,ic1,ic2,ic3,ic4]
      c[1,2,ic1,ic2,ic3,ic4] -= coeff*c[2,2,ic1,ic2,ic3,ic4]
      c[1,3,ic1,ic2,ic3,ic4] -= coeff*c[2,3,ic1,ic2,ic3,ic4]
      c[1,4,ic1,ic2,ic3,ic4] -= coeff*c[2,4,ic1,ic2,ic3,ic4]
      c[1,5,ic1,ic2,ic3,ic4] -= coeff*c[2,5,ic1,ic2,ic3,ic4]
      r[1,ir1,ir2,ir3,ir4] -= coeff*r[2,ir1,ir2,ir3,ir4]

      coeff = lhs[3,2,i]
      lhs[3,3,i] -= coeff*lhs[2,3,i]
      lhs[3,4,i] -= coeff*lhs[2,4,i]
      lhs[3,5,i] -= coeff*lhs[2,5,i]
      c[3,1,ic1,ic2,ic3,ic4] -= coeff*c[2,1,ic1,ic2,ic3,ic4]
      c[3,2,ic1,ic2,ic3,ic4] -= coeff*c[2,2,ic1,ic2,ic3,ic4]
      c[3,3,ic1,ic2,ic3,ic4] -= coeff*c[2,3,ic1,ic2,ic3,ic4]
      c[3,4,ic1,ic2,ic3,ic4] -= coeff*c[2,4,ic1,ic2,ic3,ic4]
      c[3,5,ic1,ic2,ic3,ic4] -= coeff*c[2,5,ic1,ic2,ic3,ic4]
      r[3,ir1,ir2,ir3,ir4] -= coeff*r[2,ir1,ir2,ir3,ir4]

      coeff = lhs[4,2,i]
      lhs[4,3,i] -= coeff*lhs[2,3,i]
      lhs[4,4,i] -= coeff*lhs[2,4,i]
      lhs[4,5,i] -= coeff*lhs[2,5,i]
      c[4,1,ic1,ic2,ic3,ic4] -= coeff*c[2,1,ic1,ic2,ic3,ic4]
      c[4,2,ic1,ic2,ic3,ic4] -= coeff*c[2,2,ic1,ic2,ic3,ic4]
      c[4,3,ic1,ic2,ic3,ic4] -= coeff*c[2,3,ic1,ic2,ic3,ic4]
      c[4,4,ic1,ic2,ic3,ic4] -= coeff*c[2,4,ic1,ic2,ic3,ic4]
      c[4,5,ic1,ic2,ic3,ic4] -= coeff*c[2,5,ic1,ic2,ic3,ic4]
      r[4,ir1,ir2,ir3,ir4] -= coeff*r[2,ir1,ir2,ir3,ir4]

      coeff = lhs[5,2,i]
      lhs[5,3,i] -= coeff*lhs[2,3,i]
      lhs[5,4,i] -= coeff*lhs[2,4,i]
      lhs[5,5,i] -= coeff*lhs[2,5,i]
      c[5,1,ic1,ic2,ic3,ic4] -= coeff*c[2,1,ic1,ic2,ic3,ic4]
      c[5,2,ic1,ic2,ic3,ic4] -= coeff*c[2,2,ic1,ic2,ic3,ic4]
      c[5,3,ic1,ic2,ic3,ic4] -= coeff*c[2,3,ic1,ic2,ic3,ic4]
      c[5,4,ic1,ic2,ic3,ic4] -= coeff*c[2,4,ic1,ic2,ic3,ic4]
      c[5,5,ic1,ic2,ic3,ic4] -= coeff*c[2,5,ic1,ic2,ic3,ic4]
      r[5,ir1,ir2,ir3,ir4] -= coeff*r[2,ir1,ir2,ir3,ir4]

      pivot = 1.00e0/lhs[3,3,i]
      lhs[3,4,i] = lhs[3,4,i]*pivot
      lhs[3,5,i] = lhs[3,5,i]*pivot
      c[3,1,ic1,ic2,ic3,ic4] *= pivot
      c[3,2,ic1,ic2,ic3,ic4] *= pivot
      c[3,3,ic1,ic2,ic3,ic4] *= pivot
      c[3,4,ic1,ic2,ic3,ic4] *= pivot
      c[3,5,ic1,ic2,ic3,ic4] *= pivot
      r[3,ir1,ir2,ir3,ir4] *= pivot

      coeff = lhs[1,3,i]
      lhs[1,4,i] -= coeff*lhs[3,4,i]
      lhs[1,5,i] -= coeff*lhs[3,5,i]
      c[1,1,ic1,ic2,ic3,ic4] -= coeff*c[3,1,ic1,ic2,ic3,ic4]
      c[1,2,ic1,ic2,ic3,ic4] -= coeff*c[3,2,ic1,ic2,ic3,ic4]
      c[1,3,ic1,ic2,ic3,ic4] -= coeff*c[3,3,ic1,ic2,ic3,ic4]
      c[1,4,ic1,ic2,ic3,ic4] -= coeff*c[3,4,ic1,ic2,ic3,ic4]
      c[1,5,ic1,ic2,ic3,ic4] -= coeff*c[3,5,ic1,ic2,ic3,ic4]
      r[1,ir1,ir2,ir3,ir4] -= coeff*r[3,ir1,ir2,ir3,ir4]

      coeff = lhs[2,3,i]
      lhs[2,4,i] -= coeff*lhs[3,4,i]
      lhs[2,5,i] -= coeff*lhs[3,5,i]
      c[2,1,ic1,ic2,ic3,ic4] -= coeff*c[3,1,ic1,ic2,ic3,ic4]
      c[2,2,ic1,ic2,ic3,ic4] -= coeff*c[3,2,ic1,ic2,ic3,ic4]
      c[2,3,ic1,ic2,ic3,ic4] -= coeff*c[3,3,ic1,ic2,ic3,ic4]
      c[2,4,ic1,ic2,ic3,ic4] -= coeff*c[3,4,ic1,ic2,ic3,ic4]
      c[2,5,ic1,ic2,ic3,ic4] -= coeff*c[3,5,ic1,ic2,ic3,ic4]
      r[2,ir1,ir2,ir3,ir4] -= coeff*r[3,ir1,ir2,ir3,ir4]

      coeff = lhs[4,3,i]
      lhs[4,4,i] -= coeff*lhs[3,4,i]
      lhs[4,5,i] -= coeff*lhs[3,5,i]
      c[4,1,ic1,ic2,ic3,ic4] -= coeff*c[3,1,ic1,ic2,ic3,ic4]
      c[4,2,ic1,ic2,ic3,ic4] -= coeff*c[3,2,ic1,ic2,ic3,ic4]
      c[4,3,ic1,ic2,ic3,ic4] -= coeff*c[3,3,ic1,ic2,ic3,ic4]
      c[4,4,ic1,ic2,ic3,ic4] -= coeff*c[3,4,ic1,ic2,ic3,ic4]
      c[4,5,ic1,ic2,ic3,ic4] -= coeff*c[3,5,ic1,ic2,ic3,ic4]
      r[4,ir1,ir2,ir3,ir4] -= coeff*r[3,ir1,ir2,ir3,ir4]

      coeff = lhs[5,3,i]
      lhs[5,4,i] -= coeff*lhs[3,4,i]
      lhs[5,5,i] -= coeff*lhs[3,5,i]
      c[5,1,ic1,ic2,ic3,ic4] -= coeff*c[3,1,ic1,ic2,ic3,ic4]
      c[5,2,ic1,ic2,ic3,ic4] -= coeff*c[3,2,ic1,ic2,ic3,ic4]
      c[5,3,ic1,ic2,ic3,ic4] -= coeff*c[3,3,ic1,ic2,ic3,ic4]
      c[5,4,ic1,ic2,ic3,ic4] -= coeff*c[3,4,ic1,ic2,ic3,ic4]
      c[5,5,ic1,ic2,ic3,ic4] -= coeff*c[3,5,ic1,ic2,ic3,ic4]
      r[5,ir1,ir2,ir3,ir4] -= coeff*r[3,ir1,ir2,ir3,ir4]

      pivot = 1.00e0/lhs[4,4,i]
      lhs[4,5,i] = lhs[4,5,i]*pivot
      c[4,1,ic1,ic2,ic3,ic4] *= pivot
      c[4,2,ic1,ic2,ic3,ic4] *= pivot
      c[4,3,ic1,ic2,ic3,ic4] *= pivot
      c[4,4,ic1,ic2,ic3,ic4] *= pivot
      c[4,5,ic1,ic2,ic3,ic4] *= pivot
      r[4,ir1,ir2,ir3,ir4] *= pivot

      coeff = lhs[1,4,i]
      lhs[1,5,i] -= coeff*lhs[4,5,i]
      c[1,1,ic1,ic2,ic3,ic4] -= coeff*c[4,1,ic1,ic2,ic3,ic4]
      c[1,2,ic1,ic2,ic3,ic4] -= coeff*c[4,2,ic1,ic2,ic3,ic4]
      c[1,3,ic1,ic2,ic3,ic4] -= coeff*c[4,3,ic1,ic2,ic3,ic4]
      c[1,4,ic1,ic2,ic3,ic4] -= coeff*c[4,4,ic1,ic2,ic3,ic4]
      c[1,5,ic1,ic2,ic3,ic4] -= coeff*c[4,5,ic1,ic2,ic3,ic4]
      r[1,ir1,ir2,ir3,ir4] -= coeff*r[4,ir1,ir2,ir3,ir4]

      coeff = lhs[2,4,i]
      lhs[2,5,i] -= coeff*lhs[4,5,i]
      c[2,1,ic1,ic2,ic3,ic4] -= coeff*c[4,1,ic1,ic2,ic3,ic4]
      c[2,2,ic1,ic2,ic3,ic4] -= coeff*c[4,2,ic1,ic2,ic3,ic4]
      c[2,3,ic1,ic2,ic3,ic4] -= coeff*c[4,3,ic1,ic2,ic3,ic4]
      c[2,4,ic1,ic2,ic3,ic4] -= coeff*c[4,4,ic1,ic2,ic3,ic4]
      c[2,5,ic1,ic2,ic3,ic4] -= coeff*c[4,5,ic1,ic2,ic3,ic4]
      r[2,ir1,ir2,ir3,ir4] -= coeff*r[4,ir1,ir2,ir3,ir4]

      coeff = lhs[3,4,i]
      lhs[3,5,i] -= coeff*lhs[4,5,i]
      c[3,1,ic1,ic2,ic3,ic4] -= coeff*c[4,1,ic1,ic2,ic3,ic4]
      c[3,2,ic1,ic2,ic3,ic4] -= coeff*c[4,2,ic1,ic2,ic3,ic4]
      c[3,3,ic1,ic2,ic3,ic4] -= coeff*c[4,3,ic1,ic2,ic3,ic4]
      c[3,4,ic1,ic2,ic3,ic4] -= coeff*c[4,4,ic1,ic2,ic3,ic4]
      c[3,5,ic1,ic2,ic3,ic4] -= coeff*c[4,5,ic1,ic2,ic3,ic4]
      r[3,ir1,ir2,ir3,ir4] -= coeff*r[4,ir1,ir2,ir3,ir4]

      coeff = lhs[5,4,i]
      lhs[5,5,i] -= coeff*lhs[4,5,i]
      c[5,1,ic1,ic2,ic3,ic4] -= coeff*c[4,1,ic1,ic2,ic3,ic4]
      c[5,2,ic1,ic2,ic3,ic4] -= coeff*c[4,2,ic1,ic2,ic3,ic4]
      c[5,3,ic1,ic2,ic3,ic4] -= coeff*c[4,3,ic1,ic2,ic3,ic4]
      c[5,4,ic1,ic2,ic3,ic4] -= coeff*c[4,4,ic1,ic2,ic3,ic4]
      c[5,5,ic1,ic2,ic3,ic4] -= coeff*c[4,5,ic1,ic2,ic3,ic4]
      r[5,ir1,ir2,ir3,ir4] -= coeff*r[4,ir1,ir2,ir3,ir4]


      pivot = 1.00e0/lhs[5,5,i]
      c[5,1,ic1,ic2,ic3,ic4] *= pivot
      c[5,2,ic1,ic2,ic3,ic4] *= pivot
      c[5,3,ic1,ic2,ic3,ic4] *= pivot
      c[5,4,ic1,ic2,ic3,ic4] *= pivot
      c[5,5,ic1,ic2,ic3,ic4] *= pivot
      r[5,ir1,ir2,ir3,ir4]   *= pivot

      coeff = lhs[1,5,i]
      c[1,1,ic1,ic2,ic3,ic4] -= coeff*c[5,1,ic1,ic2,ic3,ic4]
      c[1,2,ic1,ic2,ic3,ic4] -= coeff*c[5,2,ic1,ic2,ic3,ic4]
      c[1,3,ic1,ic2,ic3,ic4] -= coeff*c[5,3,ic1,ic2,ic3,ic4]
      c[1,4,ic1,ic2,ic3,ic4] -= coeff*c[5,4,ic1,ic2,ic3,ic4]
      c[1,5,ic1,ic2,ic3,ic4] -= coeff*c[5,5,ic1,ic2,ic3,ic4]
      r[1,ir1,ir2,ir3,ir4] -= coeff*r[5,ir1,ir2,ir3,ir4]

      coeff = lhs[2,5,i]
      c[2,1,ic1,ic2,ic3,ic4] -= coeff*c[5,1,ic1,ic2,ic3,ic4]
      c[2,2,ic1,ic2,ic3,ic4] -= coeff*c[5,2,ic1,ic2,ic3,ic4]
      c[2,3,ic1,ic2,ic3,ic4] -= coeff*c[5,3,ic1,ic2,ic3,ic4]
      c[2,4,ic1,ic2,ic3,ic4] -= coeff*c[5,4,ic1,ic2,ic3,ic4]
      c[2,5,ic1,ic2,ic3,ic4] -= coeff*c[5,5,ic1,ic2,ic3,ic4]
      r[2,ir1,ir2,ir3,ir4] -= coeff*r[5,ir1,ir2,ir3,ir4]

      coeff = lhs[3,5,i]
      c[3,1,ic1,ic2,ic3,ic4] -= coeff*c[5,1,ic1,ic2,ic3,ic4]
      c[3,2,ic1,ic2,ic3,ic4] -= coeff*c[5,2,ic1,ic2,ic3,ic4]
      c[3,3,ic1,ic2,ic3,ic4] -= coeff*c[5,3,ic1,ic2,ic3,ic4]
      c[3,4,ic1,ic2,ic3,ic4] -= coeff*c[5,4,ic1,ic2,ic3,ic4]
      c[3,5,ic1,ic2,ic3,ic4] -= coeff*c[5,5,ic1,ic2,ic3,ic4]
      r[3,ir1,ir2,ir3,ir4] -= coeff*r[5,ir1,ir2,ir3,ir4]

      coeff = lhs[4,5,i]
      c[4,1,ic1,ic2,ic3,ic4] -= coeff*c[5,1,ic1,ic2,ic3,ic4]
      c[4,2,ic1,ic2,ic3,ic4] -= coeff*c[5,2,ic1,ic2,ic3,ic4]
      c[4,3,ic1,ic2,ic3,ic4] -= coeff*c[5,3,ic1,ic2,ic3,ic4]
      c[4,4,ic1,ic2,ic3,ic4] -= coeff*c[5,4,ic1,ic2,ic3,ic4]
      c[4,5,ic1,ic2,ic3,ic4] -= coeff*c[5,5,ic1,ic2,ic3,ic4]
      r[4,ir1,ir2,ir3,ir4] -= coeff*r[5,ir1,ir2,ir3,ir4]

      return nothing
end



#---------------------------------------------------------------------
#---------------------------------------------------------------------

function binvrhs(lhsb, i, rhs, ir1, ir2, ir3, ir4)

      pivot = 1.00e0/lhsb[1,1,i]
      lhsb[1,2,i] *= pivot
      lhsb[1,3,i] *= pivot
      lhsb[1,4,i] *= pivot
      lhsb[1,5,i] *= pivot
      rhs[1,ir1,ir2,ir3,ir4] *= pivot

      coeff = lhsb[2,1,i]
      lhsb[2,2,i] -= coeff*lhsb[1,2,i]
      lhsb[2,3,i] -= coeff*lhsb[1,3,i]
      lhsb[2,4,i] -= coeff*lhsb[1,4,i]
      lhsb[2,5,i] -= coeff*lhsb[1,5,i]
      rhs[2,ir1,ir2,ir3,ir4] -= coeff*rhs[1,ir1,ir2,ir3,ir4]

      coeff = lhsb[3,1,i]
      lhsb[3,2,i] -= coeff*lhsb[1,2,i]
      lhsb[3,3,i] -= coeff*lhsb[1,3,i]
      lhsb[3,4,i] -= coeff*lhsb[1,4,i]
      lhsb[3,5,i] -= coeff*lhsb[1,5,i]
      rhs[3,ir1,ir2,ir3,ir4] -= coeff*rhs[1,ir1,ir2,ir3,ir4]

      coeff = lhsb[4,1,i]
      lhsb[4,2,i] -= coeff*lhsb[1,2,i]
      lhsb[4,3,i] -= coeff*lhsb[1,3,i]
      lhsb[4,4,i] -= coeff*lhsb[1,4,i]
      lhsb[4,5,i] -= coeff*lhsb[1,5,i]
      rhs[4,ir1,ir2,ir3,ir4] -= coeff*rhs[1,ir1,ir2,ir3,ir4]

      coeff = lhsb[5,1,i]
      lhsb[5,2,i] -= coeff*lhsb[1,2,i]
      lhsb[5,3,i] -= coeff*lhsb[1,3,i]
      lhsb[5,4,i] -= coeff*lhsb[1,4,i]
      lhsb[5,5,i] -= coeff*lhsb[1,5,i]
      rhs[5,ir1,ir2,ir3,ir4] -= coeff*rhs[1,ir1,ir2,ir3,ir4]


      pivot = 1.00e0/lhsb[2,2,i]
      lhsb[2,3,i] *= pivot
      lhsb[2,4,i] *= pivot
      lhsb[2,5,i] *= pivot
      rhs[2,ir1,ir2,ir3,ir4] *= pivot

      coeff = lhsb[1,2,i]
      lhsb[1,3,i] -= coeff*lhsb[2,3,i]
      lhsb[1,4,i] -= coeff*lhsb[2,4,i]
      lhsb[1,5,i] -= coeff*lhsb[2,5,i]
      rhs[1,ir1,ir2,ir3,ir4] -= coeff*rhs[2,ir1,ir2,ir3,ir4]

      coeff = lhsb[3,2,i]
      lhsb[3,3,i] -= coeff*lhsb[2,3,i]
      lhsb[3,4,i] -= coeff*lhsb[2,4,i]
      lhsb[3,5,i] -= coeff*lhsb[2,5,i]
      rhs[3,ir1,ir2,ir3,ir4] -= coeff*rhs[2,ir1,ir2,ir3,ir4]

      coeff = lhsb[4,2,i]
      lhsb[4,3,i] -= coeff*lhsb[2,3,i]
      lhsb[4,4,i] -= coeff*lhsb[2,4,i]
      lhsb[4,5,i] -= coeff*lhsb[2,5,i]
      rhs[4,ir1,ir2,ir3,ir4] -= coeff*rhs[2,ir1,ir2,ir3,ir4]

      coeff = lhsb[5,2,i]
      lhsb[5,3,i] -= coeff*lhsb[2,3,i]
      lhsb[5,4,i] -= coeff*lhsb[2,4,i]
      lhsb[5,5,i] -= coeff*lhsb[2,5,i]
      rhs[5,ir1,ir2,ir3,ir4] -= coeff*rhs[2,ir1,ir2,ir3,ir4]


      pivot = 1.00e0/lhsb[3,3,i]
      lhsb[3,4,i] *= pivot
      lhsb[3,5,i] *= pivot
      rhs[3,ir1,ir2,ir3,ir4] *= pivot

      coeff = lhsb[1,3,i]
      lhsb[1,4,i] -= coeff*lhsb[3,4,i]
      lhsb[1,5,i] -= coeff*lhsb[3,5,i]
      rhs[1,ir1,ir2,ir3,ir4] -= coeff*rhs[3,ir1,ir2,ir3,ir4]

      coeff = lhsb[2,3,i]
      lhsb[2,4,i] -= coeff*lhsb[3,4,i]
      lhsb[2,5,i] -= coeff*lhsb[3,5,i]
      rhs[2,ir1,ir2,ir3,ir4]  -= coeff*rhs[3,ir1,ir2,ir3,ir4]

      coeff = lhsb[4,3,i]
      lhsb[4,4,i] -= coeff*lhsb[3,4,i]
      lhsb[4,5,i] -= coeff*lhsb[3,5,i]
      rhs[4,ir1,ir2,ir3,ir4] -= coeff*rhs[3,ir1,ir2,ir3,ir4]

      coeff = lhsb[5,3,i]
      lhsb[5,4,i] -= coeff*lhsb[3,4,i]
      lhsb[5,5,i] -= coeff*lhsb[3,5,i]
      rhs[5,ir1,ir2,ir3,ir4] -= coeff*rhs[3,ir1,ir2,ir3,ir4]


      pivot = 1.00e0/lhsb[4,4,i]
      lhsb[4,5,i] *= pivot
      rhs[4,ir1,ir2,ir3,ir4] *= pivot

      coeff = lhsb[1,4,i]
      lhsb[1,5,i] -= coeff*lhsb[4,5,i]
      rhs[1,ir1,ir2,ir3,ir4] -= coeff*rhs[4,ir1,ir2,ir3,ir4]

      coeff = lhsb[2,4,i]
      lhsb[2,5,i] -= coeff*lhsb[4,5,i]
      rhs[2,ir1,ir2,ir3,ir4] -= coeff*rhs[4,ir1,ir2,ir3,ir4]

      coeff = lhsb[3,4,i]
      lhsb[3,5,i] -= coeff*lhsb[4,5,i]
      rhs[3,ir1,ir2,ir3,ir4] -= coeff*rhs[4,ir1,ir2,ir3,ir4]

      coeff = lhsb[5,4,i]
      lhsb[5,5,i] -= coeff*lhsb[4,5,i]
      rhs[5,ir1,ir2,ir3,ir4] -= coeff*rhs[4,ir1,ir2,ir3,ir4]
      

      pivot = 1.00e0/lhsb[5,5,i]
      rhs[5,ir1,ir2,ir3,ir4] *= pivot

      coeff = lhsb[1,5,i]
      rhs[1,ir1,ir2,ir3,ir4] -= coeff*rhs[5,ir1,ir2,ir3,ir4]

      coeff = lhsb[2,5,i]
      rhs[2,ir1,ir2,ir3,ir4] -= coeff*rhs[5,ir1,ir2,ir3,ir4]

      coeff = lhsb[3,5,i]
      rhs[3,ir1,ir2,ir3,ir4] -= coeff*rhs[5,ir1,ir2,ir3,ir4]

      coeff = lhsb[4,5,i]
      rhs[4,ir1,ir2,ir3,ir4] -= coeff*rhs[5,ir1,ir2,ir3,ir4]

      return nothing
end



