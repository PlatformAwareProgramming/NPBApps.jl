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


#---------------------------------------------------------------------
#  verification routine                         
#---------------------------------------------------------------------

 function verify(xcr, xce, xci, class, dt)

         xcrref = Array{Float64}(undef, 5)
         xceref = Array{Float64}(undef, 5)
         xcrdif = Array{Float64}(undef, 5)
         xcedif = Array{Float64}(undef, 5)

#---------------------------------------------------------------------
#   tolerance level
#---------------------------------------------------------------------
        epsilon = 1.0e-08

        verified = true

        for m = 1:5
           xcrref[m] = 1.0
           xceref[m] = 1.0
        end
        xciref = 1.0

        if class == CLASS_S

         dtref = 5.0e-1
#---------------------------------------------------------------------
#   Reference values of RMS-norms of residual, for the (12X12X12) grid,
#   after 50 time steps, with  DT = 5.0d-01
#---------------------------------------------------------------------
         xcrref[1] = 1.6196343210976702e-02
         xcrref[2] = 2.1976745164821318e-03
         xcrref[3] = 1.5179927653399185e-03
         xcrref[4] = 1.5029584435994323e-03
         xcrref[5] = 3.4264073155896461e-02

#---------------------------------------------------------------------
#   Reference values of RMS-norms of solution error, for the (12X12X12) grid,
#   after 50 time steps, with  DT = 5.0d-01
#---------------------------------------------------------------------
         xceref[1] = 6.4223319957960924e-04
         xceref[2] = 8.4144342047347926e-05
         xceref[3] = 5.8588269616485186e-05
         xceref[4] = 5.8474222595157350e-05
         xceref[5] = 1.3103347914111294e-03

#---------------------------------------------------------------------
#   Reference value of surface integral, for the (12X12X12) grid,
#   after 50 time steps, with DT = 5.0d-01
#---------------------------------------------------------------------
         xciref = 7.8418928865937083e+00

        elseif class == CLASS_W

           dtref = 1.5e-3
#---------------------------------------------------------------------
#   Reference values of RMS-norms of residual, for the (33x33x33) grid,
#   after 300 time steps, with  DT = 1.5d-3
#---------------------------------------------------------------------
           xcrref[1] =   0.1236511638192e+02
           xcrref[2] =   0.1317228477799e+01
           xcrref[3] =   0.2550120713095e+01
           xcrref[4] =   0.2326187750252e+01
           xcrref[5] =   0.2826799444189e+02


#---------------------------------------------------------------------
#   Reference values of RMS-norms of solution error, for the (33X33X33) grid,
#---------------------------------------------------------------------
           xceref[1] =   0.4867877144216e+00
           xceref[2] =   0.5064652880982e-01
           xceref[3] =   0.9281818101960e-01
           xceref[4] =   0.8570126542733e-01
           xceref[5] =   0.1084277417792e+01


#---------------------------------------------------------------------
#   Reference value of surface integral, for the (33X33X33) grid,
#   after 300 time steps, with  DT = 1.5d-3
#---------------------------------------------------------------------
           xciref    =   0.1161399311023e+02

        elseif class == CLASS_A

           dtref = 2.0e+0
#---------------------------------------------------------------------
#   Reference values of RMS-norms of residual, for the (64X64X64) grid,
#   after 250 time steps, with  DT = 2.0d+00
#---------------------------------------------------------------------
         xcrref[1] = 7.7902107606689367e+02
         xcrref[2] = 6.3402765259692870e+01
         xcrref[3] = 1.9499249727292479e+02
         xcrref[4] = 1.7845301160418537e+02
         xcrref[5] = 1.8384760349464247e+03

#---------------------------------------------------------------------
#   Reference values of RMS-norms of solution error, for the (64X64X64) grid,
#   after 250 time steps, with  DT = 2.0d+00
#---------------------------------------------------------------------
         xceref[1] = 2.9964085685471943e+01
         xceref[2] = 2.8194576365003349e+00
         xceref[3] = 7.3473412698774742e+00
         xceref[4] = 6.7139225687777051e+00
         xceref[5] = 7.0715315688392578e+01

#---------------------------------------------------------------------
#   Reference value of surface integral, for the (64X64X64) grid,
#   after 250 time steps, with DT = 2.0d+00
#---------------------------------------------------------------------
         xciref = 2.6030925604886277e+01

        elseif class == CLASS_B

           dtref = 2.0e+0

#---------------------------------------------------------------------
#   Reference values of RMS-norms of residual, for the (102X102X102) grid,
#   after 250 time steps, with  DT = 2.0d+00
#---------------------------------------------------------------------
         xcrref[1] = 3.5532672969982736e+03
         xcrref[2] = 2.6214750795310692e+02
         xcrref[3] = 8.8333721850952190e+02
         xcrref[4] = 7.7812774739425265e+02
         xcrref[5] = 7.3087969592545314e+03

#---------------------------------------------------------------------
#   Reference values of RMS-norms of solution error, for the (102X102X102) 
#   grid, after 250 time steps, with  DT = 2.0d+00
#---------------------------------------------------------------------
         xceref[1] = 1.1401176380212709e+02
         xceref[2] = 8.1098963655421574e+00
         xceref[3] = 2.8480597317698308e+01
         xceref[4] = 2.5905394567832939e+01
         xceref[5] = 2.6054907504857413e+02

#---------------------------------------------------------------------
#   Reference value of surface integral, for the (102X102X102) grid,
#   after 250 time steps, with DT = 2.0d+00
#---------------------------------------------------------------------
         xciref = 4.7887162703308227e+01

        elseif class == CLASS_C

           dtref = 2.0e+0

#---------------------------------------------------------------------
#   Reference values of RMS-norms of residual, for the (162X162X162) grid,
#   after 250 time steps, with  DT = 2.0d+00
#---------------------------------------------------------------------
         xcrref[1] = 1.03766980323537846e+04
         xcrref[2] = 8.92212458801008552e+02
         xcrref[3] = 2.56238814582660871e+03
         xcrref[4] = 2.19194343857831427e+03
         xcrref[5] = 1.78078057261061185e+04

#---------------------------------------------------------------------
#   Reference values of RMS-norms of solution error, for the (162X162X162) 
#   grid, after 250 time steps, with  DT = 2.0d+00
#---------------------------------------------------------------------
         xceref[1] = 2.15986399716949279e+02
         xceref[2] = 1.55789559239863600e+01
         xceref[3] = 5.41318863077207766e+01
         xceref[4] = 4.82262643154045421e+01
         xceref[5] = 4.55902910043250358e+02

#---------------------------------------------------------------------
#   Reference value of surface integral, for the (162X162X162) grid,
#   after 250 time steps, with DT = 2.0d+00
#---------------------------------------------------------------------
         xciref = 6.66404553572181300e+01

        elseif class == CLASS_D

           dtref = 1.0e+0

#---------------------------------------------------------------------
#   Reference values of RMS-norms of residual, for the (408X408X408) grid,
#   after 300 time steps, with  DT = 1.0d+00
#---------------------------------------------------------------------
         xcrref[1] = 0.4868417937025e+05
         xcrref[2] = 0.4696371050071e+04
         xcrref[3] = 0.1218114549776e+05
         xcrref[4] = 0.1033801493461e+05
         xcrref[5] = 0.7142398413817e+05

#---------------------------------------------------------------------
#   Reference values of RMS-norms of solution error, for the (408X408X408) 
#   grid, after 300 time steps, with  DT = 1.0d+00
#---------------------------------------------------------------------
         xceref[1] = 0.3752393004482e+03
         xceref[2] = 0.3084128893659e+02
         xceref[3] = 0.9434276905469e+02
         xceref[4] = 0.8230686681928e+02
         xceref[5] = 0.7002620636210e+03

#---------------------------------------------------------------------
#   Reference value of surface integral, for the (408X408X408) grid,
#   after 300 time steps, with DT = 1.0d+00
#---------------------------------------------------------------------
         xciref =    0.8334101392503e+02

        elseif class == CLASS_E

           dtref = 0.5e+0

#---------------------------------------------------------------------
#   Reference values of RMS-norms of residual, for the (1020X1020X1020) grid,
#   after 300 time steps, with  DT = 0.5d+00
#---------------------------------------------------------------------
         xcrref[1] = 0.2099641687874e+06
         xcrref[2] = 0.2130403143165e+05
         xcrref[3] = 0.5319228789371e+05
         xcrref[4] = 0.4509761639833e+05
         xcrref[5] = 0.2932360006590e+06

#---------------------------------------------------------------------
#   Reference values of RMS-norms of solution error, for the (1020X1020X1020) 
#   grid, after 300 time steps, with  DT = 0.5d+00
#---------------------------------------------------------------------
         xceref[1] = 0.4800572578333e+03
         xceref[2] = 0.4221993400184e+02
         xceref[3] = 0.1210851906824e+03
         xceref[4] = 0.1047888986770e+03
         xceref[5] = 0.8363028257389e+03

#---------------------------------------------------------------------
#   Reference value of surface integral, for the (1020X1020X1020) grid,
#   after 300 time steps, with DT = 0.5d+00
#---------------------------------------------------------------------
         xciref =    0.9512163272273e+02

        elseif class == CLASS_F

           dtref = 0.2e+0

#---------------------------------------------------------------------
#   Reference values of RMS-norms of residual, for the (2560X2560X2560) grid,
#   after 300 time steps, with  DT = 0.2d+00
#---------------------------------------------------------------------
         xcrref[1] = 0.8505125358152e+06
         xcrref[2] = 0.8774655318044e+05
         xcrref[3] = 0.2167258198851e+06
         xcrref[4] = 0.1838245257371e+06
         xcrref[5] = 0.1175556512415e+07

#---------------------------------------------------------------------
#   Reference values of RMS-norms of solution error, for the (2560X2560X2560)
#   grid, after 300 time steps, with  DT = 0.2d+00
#---------------------------------------------------------------------
         xceref[1] = 0.5293914132486e+03
         xceref[2] = 0.4784861621068e+02
         xceref[3] = 0.1337701281659e+03
         xceref[4] = 0.1154215049655e+03
         xceref[5] = 0.8956266851467e+03

#---------------------------------------------------------------------
#   Reference value of surface integral, for the (2560X2560X2560) grid,
#   after 300 time steps, with DT = 0.2d+00
#---------------------------------------------------------------------
         xciref =    0.1002509436546e+03

        else

           verified = false

        end

#---------------------------------------------------------------------
#    verification test for residuals if gridsize is one of 
#    the defined grid sizes above (class .ne. 'U')
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#    Compute the difference of solution values and the known reference values.
#---------------------------------------------------------------------
        for m = 1:5
           xcrdif[m] = abs((xcr[m]-xcrref[m])/xcrref[m])
           xcedif[m] = abs((xce[m]-xceref[m])/xceref[m])
                
        end
        xcidif = abs((xci - xciref)/xciref)


#---------------------------------------------------------------------
#    Output the comparison of computed results to known cases.
#---------------------------------------------------------------------

        if class != CLASS_UNDEFINED
           @printf(stdout, "\n Verification being performed for class %s\n", class)
           @printf(stdout, " Accuracy setting for epsilon = %20.13E\n", epsilon)
           verified = (abs(dt-dtref) <= epsilon)
           if !verified
              class = CLASS_UNDEFINED
              @printf(stdout, " DT does not match the reference value of %15.8E\n", dtref)
           end
        else
           @printf(stdout, " Unknown class\n", )
        end


        if class != CLASS_UNDEFINED
           @printf(stdout, " Comparison of RMS-norms of residual\n", )
        else
           @printf(stdout, " RMS-norms of residual\n", )
        end

        for m = 1:5
           if class == CLASS_UNDEFINED
              @printf(stdout, "          %2i  %20.13E\n", m, xcr[m])
           elseif xcrdif[m] != NaN && xcrdif[m] <= epsilon
              @printf(stdout, "          %2i  %20.13E%20.13E%20.13E\n", m, xcr[m], xcrref[m], xcrdif[m])
           else
              verified = false
              @printf(stdout, " FAILURE: %2i  %20.13E%20.13E%20.13E\n", m, xcr[m], xcrref[m], xcrdif[m])
           end
        end

        if class != CLASS_UNDEFINED
           @printf(stdout, " Comparison of RMS-norms of solution error\n", )
        else
           @printf(stdout, " RMS-norms of solution error\n", )
        end

        for m = 1:5
           if class == CLASS_UNDEFINED
              @printf(stdout, "          %2i  %20.13E\n", m, xce[m])
           elseif xcedif[m] != NaN && xcedif[m] <= epsilon
              @printf(stdout, "          %2i  %20.13E%20.13E%20.13E\n", m, xce[m], xceref[m], xcedif[m])
           else
              verified = false
              @printf(stdout, " FAILURE: %2i  %20.13E%20.13E%20.13E\n", m, xce[m], xceref[m], xcedif[m])
           end
        end

        if class != CLASS_UNDEFINED
           @printf(stdout, " Comparison of surface integral\n", )
        else
           @printf(stdout, " Surface integral\n", )
        end

        if class == CLASS_UNDEFINED
           @printf(stdout, "              %20.13E\n", xci)
        elseif xcidif <= epsilon
           @printf(stdout, "              %20.13E%20.13E%20.13E\n", xci, xciref, xcidif)
        else
           verified = false
           @printf(stdout, " FAILURE:     %20.13E%20.13E%20.13E\n", xci, xciref, xcidif)
        end

        if class == CLASS_UNDEFINED
           @printf(stdout, " No reference values provided\n", )
           @printf(stdout, " No verification performed\n", )
        elseif verified
           @printf(stdout, " Verification Successful\n", )
        else
           @printf(stdout, " Verification failed\n", )
        end

        return verified

end
