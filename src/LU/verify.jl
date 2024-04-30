#---------------------------------------------------------------------
#  set problem class based on problem size
#---------------------------------------------------------------------

 function set_class(itmax, nx0, ny0, nz0)

        if (nx0 == 12) && (ny0 == 12) && (nz0 == 12) && (itmax == 50)
           class = CLASS_S
        elseif (nx0 == 33) && (ny0 == 33) && (nz0 == 33) && (itmax  == 300)
           class = CLASS_W   #SPEC95fp size
        elseif (nx0 == 64) &&(ny0 == 64) &&(nz0 == 64) && (itmax  == 250)
           class = CLASS_A
        elseif (nx0 == 102) && (ny0 == 102) && (nz0 == 102) && (itmax  == 250)
           class = CLASS_B
        elseif (nx0 == 162) && (ny0 == 162) && (nz0 == 162) && (itmax  == 250)
           class = CLASS_C
        elseif (nx0 == 408) &&(ny0 == 408) && (nz0 == 408) && (itmax  == 300)
           class = CLASS_D
        elseif (nx0 == 1020) && (ny0 == 1020) && (nz0 == 1020) && (itmax  == 300)
           class = CLASS_E
        elseif (nx0 == 2560) && (ny0 == 2560) && (nz0 == 2560) && (itmax  == 300)
           class = CLASS_F
        else
           class = CLASS_UNDEFINED
        end

        return class
end

# AT CLUSTER MASTER
function putVerificationReport(xcr_node, xce_node, xci_node)
   remotecall(putVerificationReport, 1, clusterid, xcr_node, xce_node, xci_node; role=:worker)
end

# AT DRIVER MASTER
function putVerificationReport(clusterid, xcr_cluster, xce_cluster, xci_cluster)
   lock(verified_signal[clusterid+1])
   try
      lock(lk_update_verify)
      try
         xcr .= xcr + xcr_cluster
         xce .= xce + xce_cluster
         xci[] = xci[] + xci_cluster
      finally
         unlock(lk_update_verify)
      end
      verified_signal_condition[clusterid+1] = true
      notify(verified_signal[clusterid+1])
   finally
      unlock(verified_signal[clusterid+1])
   end
end


#---------------------------------------------------------------------
#  verification routine                         
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#---------------------------------------------------------------------
# AT DRIVER MASTER
function verify(class::CLASS, dt, itmax)

         global verified_signal_condition = Array{Bool}(undef, num_clusters)
         global verified_signal = Array{Threads.Condition}(undef, num_clusters)
         for i in 1:num_clusters
            verified_signal_condition[i] = false
            verified_signal[i] = Threads.Condition()             
         end

         global lk_update_verify = ReentrantLock()

         xcrref = Array{FloatType}(undef, 5)
         xceref = Array{FloatType}(undef, 5)
         xcrdif = Array{FloatType}(undef, 5)
         xcedif = Array{FloatType}(undef, 5)

#---------------------------------------------------------------------
#   tolerance level
#---------------------------------------------------------------------
        epsilon = FloatType == Float64 ? 1.0e-08 : 1.0e-04

        for m = 1:5
           xcrref[m] = 1.0
           xceref[m] = 1.0
        end
        xciref = 1.0

        if class == CLASS_S

           dtref = 5.0e-1
           itmaxref = 50
#---------------------------------------------------------------------
#   Reference values of RMS-norms of residual
#---------------------------------------------------------------------
           xcrref[1] = 0.3778579699366e+01
           xcrref[2] = 0.3120418698065e+00
           xcrref[3] = 0.8386213407018e+00
           xcrref[4] = 0.4452165980488e+00
           xcrref[5] = 0.7808656756434e+01

#---------------------------------------------------------------------
#   Reference values of RMS-norms of solution error
#---------------------------------------------------------------------
           xceref[1] = 0.2429480066305e+02
           xceref[2] = 0.9072817470024e+01
           xceref[3] = 0.1032621825644e+02
           xceref[4] = 0.9256791727838e+01
           xceref[5] = 0.1639045777714e+02

#---------------------------------------------------------------------
#   Reference value of surface integral
#---------------------------------------------------------------------
           xciref    = 0.4964435445706e+02

        elseif class == CLASS_W

           dtref = 1.5e-3
           itmaxref = 300

#---------------------------------------------------------------------
#   Reference values of RMS-norms of residual
#---------------------------------------------------------------------
           xcrref[1] = 0.8285060230339e+03
           xcrref[2] = 0.5753415004693e+02
           xcrref[3] = 0.2023477570531e+03
           xcrref[4] = 0.1586275182502e+03
           xcrref[5] = 0.1733925947816e+04

#---------------------------------------------------------------------
#   Reference values of RMS-norms of solution error
#---------------------------------------------------------------------
           xceref[1] = 0.7514670702651e+02
           xceref[2] = 0.9776687033238e+01
           xceref[3] = 0.2141754291209e+02
           xceref[4] = 0.1685405918675e+02
           xceref[5] = 0.1856944519722e+03

#---------------------------------------------------------------------
#   Reference value of surface integral
#---------------------------------------------------------------------
           xciref    = 0.3781055348911e+03

        elseif class == CLASS_A

           dtref = 2.0e+0
           itmaxref = 250

#---------------------------------------------------------------------
#   Reference values of RMS-norms of residual
#---------------------------------------------------------------------
           xcrref[1] = 0.1131574877175e+04
           xcrref[2] = 0.7965206944742e+02
           xcrref[3] = 0.2705587159526e+03
           xcrref[4] = 0.2129567530746e+03
           xcrref[5] = 0.2260584655432e+04

#---------------------------------------------------------------------
#   Reference values of RMS-norms of solution error
#---------------------------------------------------------------------
           xceref[1] = 0.1115694885382e+03
           xceref[2] = 0.1089257673798e+02
           xceref[3] = 0.2905379922066e+02
           xceref[4] = 0.2216126755530e+02
           xceref[5] = 0.2501762341026e+03

#---------------------------------------------------------------------
#   Reference value of surface integral
#---------------------------------------------------------------------
           xciref    = 0.5904992211511e+03

        elseif class == CLASS_B

           dtref = 2.0e+0
           itmaxref = 250

#---------------------------------------------------------------------
#   Reference values of RMS-norms of residual
#---------------------------------------------------------------------
           xcrref[1] = 0.1734656959567e+05
           xcrref[2] = 0.1238977748533e+04
           xcrref[3] = 0.4123885357100e+04
           xcrref[4] = 0.3613705834056e+04
           xcrref[5] = 0.3531187871586e+05

#---------------------------------------------------------------------
#   Reference values of RMS-norms of solution error
#---------------------------------------------------------------------
           xceref[1] = 0.1781612313296e+04
           xceref[2] = 0.1177971120769e+03
           xceref[3] = 0.4233792871440e+03
           xceref[4] = 0.3577260438230e+03
           xceref[5] = 0.3659958544012e+04

#---------------------------------------------------------------------
#   Reference value of surface integral
#---------------------------------------------------------------------
           xciref    = 0.6107041476456e+04

        elseif class == CLASS_C

           dtref = 2.0e+0
           itmaxref = 250

#---------------------------------------------------------------------
#   Reference values of RMS-norms of residual
#---------------------------------------------------------------------
           xcrref[1] = 0.4108743427233e+05
           xcrref[2] = 0.3439004802235e+04
           xcrref[3] = 0.9961331392486e+04
           xcrref[4] = 0.8321426758084e+04
           xcrref[5] = 0.7463792419218e+05

#---------------------------------------------------------------------
#   Reference values of RMS-norms of solution error
#---------------------------------------------------------------------
           xceref[1] = 0.3429276307955e+04
           xceref[2] = 0.2336680861825e+03
           xceref[3] = 0.8216363109621e+03
           xceref[4] = 0.7143809828225e+03
           xceref[5] = 0.7057470798773e+04

#---------------------------------------------------------------------
#   Reference value of surface integral
#---------------------------------------------------------------------
            xciref    = 0.1125826349653e+05

        elseif class == CLASS_D

           dtref = 1.0e+0
           itmaxref = 300

#---------------------------------------------------------------------
#   Reference values of RMS-norms of residual
#---------------------------------------------------------------------
           xcrref[1] = 0.3282253166388e+06
           xcrref[2] = 0.3490781637713e+05
           xcrref[3] = 0.8610311978292e+05
           xcrref[4] = 0.7004896022603e+05
           xcrref[5] = 0.4546838584391e+06

#---------------------------------------------------------------------
#   Reference values of RMS-norms of solution error
#---------------------------------------------------------------------
           xceref[1] = 0.6620775619126e+04
           xceref[2] = 0.5229798207352e+03
           xceref[3] = 0.1620218261697e+04
           xceref[4] = 0.1404783445006e+04
           xceref[5] = 0.1222629805121e+05

#---------------------------------------------------------------------
#   Reference value of surface integral
#---------------------------------------------------------------------
           xciref    = 0.2059421629621e+05

        elseif class == CLASS_E

           dtref = 0.5e+0
           itmaxref = 300

#---------------------------------------------------------------------
#   Reference values of RMS-norms of residual
#---------------------------------------------------------------------
           xcrref[1] = 0.1539988626779e+07
           xcrref[2] = 0.1742224758490e+06
           xcrref[3] = 0.4153598861059e+06
           xcrref[4] = 0.3468381400447e+06
           xcrref[5] = 0.2054406022038e+07

#---------------------------------------------------------------------
#   Reference values of RMS-norms of solution error
#---------------------------------------------------------------------
           xceref[1] = 0.8021145134635e+04
           xceref[2] = 0.6932079823791e+03
           xceref[3] = 0.1998959591111e+04
           xceref[4] = 0.1725962639357e+04
           xceref[5] = 0.1389447024442e+05

#---------------------------------------------------------------------
#   Reference value of surface integral
#---------------------------------------------------------------------
           xciref    = 0.2334131124791e+05

        elseif class == CLASS_F

         dtref = 0.2e+0
         itmaxref = 300

#---------------------------------------------------------------------
#   Reference values of RMS-norms of residual
#---------------------------------------------------------------------
         xcrref[1] = 0.7116424271317e+07
         xcrref[2] = 0.8159357680842e+06
         xcrref[3] = 0.1930561069782e+07
         xcrref[4] = 0.1633447037519e+07
         xcrref[5] = 0.9417323380798e+07

#---------------------------------------------------------------------
#   Reference values of RMS-norms of solution error
#---------------------------------------------------------------------
         xceref[1] = 0.8648720989200e+04
         xceref[2] = 0.7774221260694e+03
         xceref[3] = 0.2175462599498e+04
         xceref[4] = 0.1875280641999e+04
         xceref[5] = 0.1457903413233e+05

#---------------------------------------------------------------------
#   Reference value of surface integral
#---------------------------------------------------------------------
         xciref    = 0.2448986519022e+05

         if itmax == 30

            itmaxref = 30
            xcrref[1] = 0.3814950058736e+08
            xcrref[2] = 0.4280439009977e+07
            xcrref[3] = 0.1016353864923e+08
            xcrref[4] = 0.8627208852987e+07
            xcrref[5] = 0.5024448179760e+08
 
            xceref[1] = 0.8903253221139e+04
            xceref[2] = 0.8129462858441e+03
            xceref[3] = 0.2248648703838e+04
            xceref[4] = 0.1937258920446e+04
            xceref[5] = 0.1485251162647e+05
 
            xciref    = 0.2792087395236e+05
 
         end
      else

            dtref = 0.0e+0
            itmaxref = 0
            verified = false
 
      end

#---------------------------------------------------------------------
#    verification test for residuals if gridsize is one of 
#    the defined grid sizes above (class .ne. 'U')
#---------------------------------------------------------------------

        for i = 1:num_clusters
            lock(verified_signal[i])
            try
               while !verified_signal_condition[i]
                  wait(verified_signal[i])
               end
            finally
               unlock(verified_signal[i])
            end
        end
     
        verified = true

#---------------------------------------------------------------------
#    Compute the difference of solution values and the known reference values.
#---------------------------------------------------------------------
        xcrdif = abs.((xcr .- xcrref)./xcrref)
        xcedif = abs.((xce .- xceref)./xceref)
        xcidif = abs((xci[] - xciref)/xciref)

#---------------------------------------------------------------------
#    Output the comparison of computed results to known cases.
#---------------------------------------------------------------------

         @printf(stdout, "\n Verification being performed for class %s\n", class)
         @printf(stdout, " Accuracy setting for epsilon = %20.13E\n", epsilon)
         verified = (abs(dt-dtref) <= epsilon)
         if !verified
            class = CLASS_UNDEFINED
            @printf(stdout, " DT does not match the reference value of %15.8E\n", dtref)
         elseif itmax != itmaxref
            verified = false
            @printf(stdout, " ITMAX does not match the reference value of %5i\n", itmaxref)
         end

         @printf(stdout, " Comparison of RMS-norms of residual\n", )

        for m = 1:5
           if xcrdif[m] != NaN && xcrdif[m] <= epsilon
              @printf(stdout, "          %2i  %20.13E%20.13E%20.13E\n", m, xcr[m], xcrref[m], xcrdif[m])
           else
              verified = false
              @printf(stdout, " FAILURE: %2i  %20.13E%20.13E%20.13E\n", m, xcr[m], xcrref[m], xcrdif[m])
           end
        end

        @printf(stdout, " Comparison of RMS-norms of solution error\n", )

        for m = 1:5
           if xcedif[m] != NaN && xcedif[m] <= epsilon
              @printf(stdout, "          %2i  %20.13E%20.13E%20.13E\n", m, xce[m], xceref[m], xcedif[m])
           else
              verified = false
              @printf(stdout, " FAILURE: %2i  %20.13E%20.13E%20.13E\n", m, xce[m], xceref[m], xcedif[m])
           end
        end

        @printf(stdout, " Comparison of surface integral\n", )

        if xcidif <= epsilon
           @printf(stdout, "              %20.13E%20.13E%20.13E\n", xci[], xciref, xcidif)
        else
           verified = false
           @printf(stdout, " FAILURE:     %20.13E%20.13E%20.13E\n", xci[], xciref, xcidif)
        end

        if verified
           @printf(stdout, " Verification Successful\n", )
        else
           @printf(stdout, " Verification failed\n", )
        end

        return verified
end


function verify(xcr, xce, xci)

   if node == root
      remotecall(putVerificationReport, 1, xcr, xce, xci; role = :worker)
   end

end