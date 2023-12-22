#---------------------------------------------------------------------
#
#   compute the exact solution at (i,j,k)
#
#---------------------------------------------------------------------

function exact(i, j, k, u000ijk)

      xi  = ( float( i - 1 ) ) / ( nx0 - 1 )
      eta  = ( float( j - 1 ) ) / ( ny0 - 1 )
      zeta = ( float( k - 1 ) ) / ( nz - 1 )

      for m = 1:5
         u000ijk[m] =  ce[m, 1] +
                       ce[m, 2] * xi+
                       ce[m, 3] * eta+
                       ce[m, 4] * zeta+
                       ce[m, 5] * xi * xi+
                       ce[m, 6] * eta * eta+
                       ce[m, 7] * zeta * zeta+
                       ce[m, 8] * xi * xi * xi+
                       ce[m, 9] * eta * eta * eta+
                       ce[m, 10] * zeta * zeta * zeta+
                       ce[m, 11] * xi * xi * xi * xi+
                       ce[m, 12] * eta * eta * eta * eta+
                       ce[m, 13] * zeta * zeta * zeta * zeta
      end

      return nothing
end
