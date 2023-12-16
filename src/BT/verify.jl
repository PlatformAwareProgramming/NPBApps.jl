using FortranFiles
using OffsetArrays
using Parameters
using Printf


#---------------------------------------------------------------------
#---------------------------------------------------------------------

        function set_class(no_time_steps, class)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#  set problem class based on problem size
#---------------------------------------------------------------------

#        use bt_data
#        implicit none

#        integer no_time_steps
#        character class


        if (predecessor[1]  == 12     ) &&(
             predecessor[2]  == 12     ) &&(
             predecessor[3]  == 12     ) &&(
             no_time_steps   == 60    )

           class = "S"

        elseif (predecessor[1] == 24) &&(
                 predecessor[2] == 24) &&(
                 predecessor[3] == 24) &&(
                 no_time_steps  == 200)

           class = "W"

        elseif (predecessor[1] == 64) &&(
                 predecessor[2] == 64) &&(
                 predecessor[3] == 64) &&(
                 no_time_steps  == 200)

           class = "A"

        elseif (predecessor[1] == 102) &&(
                 predecessor[2] == 102) &&(
                 predecessor[3] == 102) &&(
                 no_time_steps  == 200)

           class = "B"

        elseif (predecessor[1] == 162) &&(
                 predecessor[2] == 162) &&(
                 predecessor[3] == 162) &&(
                 no_time_steps  == 200)

           class = "C"

        elseif (predecessor[1] == 408) &&(
                 predecessor[2] == 408) &&(
                 predecessor[3] == 408) &&(
                 no_time_steps  == 250)

           class = "D"

        elseif (predecessor[1] == 1020) &&(
                 predecessor[2] == 1020) &&(
                 predecessor[3] == 1020) &&(
                 no_time_steps  == 250)

           class = "E"

        elseif (predecessor[1] == 2560) &&(
                 predecessor[2] == 2560) &&(
                 predecessor[3] == 2560) &&(
                 no_time_steps  == 250)

           class = "F"

        else

           class = "U"

        end

        return nothing
        end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

        function verify(class, verified)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#  verification routine                         
#---------------------------------------------------------------------


#        use bt_data
#        use mpinpb

#        implicit none

#        DOUBLEPRECISION xcrref[5],xceref[5],xcrdif[5],xcedif[5],  
#                         epsilon, xce[5], xcr[5], dtref
#        integer m
#        character class
#        logical verified

#---------------------------------------------------------------------
#   tolerance level
#---------------------------------------------------------------------
        epsilon = 1.0e-08
        verified = true

#---------------------------------------------------------------------
#   compute the error norm and the residual norm, and exit if not printing
#---------------------------------------------------------------------

        if iotype != 0
           timer_start(t_iov)
           accumulate_norms(xce)
           timer_stop(t_iov)
        else
           error_norm(xce)
        end

        copy_faces

        rhs_norm(xcr)

        for m = 1:5
           xcr[m] = xcr[m] / dt
        end

        if (node != 0) return end

        for m = 1:5
           xcrref[m] = 1.0
           xceref[m] = 1.0
        end

#---------------------------------------------------------------------
#    reference data for 12X12X12 grids after 60 time steps, with DT = 1.0d-02
#---------------------------------------------------------------------
        if class == "S"

           dtref = 1.0e-2

#---------------------------------------------------------------------
#  Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
         xcrref[1] = 1.7034283709541311e-01
         xcrref[2] = 1.2975252070034097e-02
         xcrref[3] = 3.2527926989486055e-02
         xcrref[4] = 2.6436421275166801e-02
         xcrref[5] = 1.9211784131744430e-01

#---------------------------------------------------------------------
#  Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------

         if iotype == 0
           xceref[1] = 4.9976913345811579e-04
           xceref[2] = 4.5195666782961927e-05
           xceref[3] = 7.3973765172921357e-05
           xceref[4] = 7.3821238632439731e-05
           xceref[5] = 8.9269630987491446e-04
         else
           xceref[1] = 0.1149036328945e+02
           xceref[2] = 0.9156788904727e+00
           xceref[3] = 0.2857899428614e+01
           xceref[4] = 0.2598273346734e+01
           xceref[5] = 0.2652795397547e+02
         end

#---------------------------------------------------------------------
#    reference data for 24X24X24 grids after 200 time steps, with DT = 0.8d-3
#---------------------------------------------------------------------
        elseif class == "W"

           dtref = 0.8e-3

#---------------------------------------------------------------------
#  Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.1125590409344e+03
           xcrref[2] = 0.1180007595731e+02
           xcrref[3] = 0.2710329767846e+02
           xcrref[4] = 0.2469174937669e+02
           xcrref[5] = 0.2638427874317e+03

#---------------------------------------------------------------------
#  Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------

         if iotype == 0
           xceref[1] = 0.4419655736008e+01
           xceref[2] = 0.4638531260002e+00
           xceref[3] = 0.1011551749967e+01
           xceref[4] = 0.9235878729944e+00
           xceref[5] = 0.1018045837718e+02
         else
           xceref[1] = 0.6729594398612e+02
           xceref[2] = 0.5264523081690e+01
           xceref[3] = 0.1677107142637e+02
           xceref[4] = 0.1508721463436e+02
           xceref[5] = 0.1477018363393e+03
         end


#---------------------------------------------------------------------
#    reference data for 64X64X64 grids after 200 time steps, with DT = 0.8d-3
#---------------------------------------------------------------------
        elseif class == "A"

           dtref = 0.8e-3

#---------------------------------------------------------------------
#  Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
         xcrref[1] = 1.0806346714637264e+02
         xcrref[2] = 1.1319730901220813e+01
         xcrref[3] = 2.5974354511582465e+01
         xcrref[4] = 2.3665622544678910e+01
         xcrref[5] = 2.5278963211748344e+02

#---------------------------------------------------------------------
#  Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------

         if iotype == 0
           xceref[1] = 4.2348416040525025e+00
           xceref[2] = 4.4390282496995698e-01
           xceref[3] = 9.6692480136345650e-01
           xceref[4] = 8.8302063039765474e-01
           xceref[5] = 9.7379901770829278e+00
         else
           xceref[1] = 0.6482218724961e+02
           xceref[2] = 0.5066461714527e+01
           xceref[3] = 0.1613931961359e+02
           xceref[4] = 0.1452010201481e+02
           xceref[5] = 0.1420099377681e+03
         end

#---------------------------------------------------------------------
#    reference data for 102X102X102 grids after 200 time steps,
#    with DT = 3.0d-04
#---------------------------------------------------------------------
        elseif class == "B"

           dtref = 3.0e-4

#---------------------------------------------------------------------
#  Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
         xcrref[1] = 1.4233597229287254e+03
         xcrref[2] = 9.9330522590150238e+01
         xcrref[3] = 3.5646025644535285e+02
         xcrref[4] = 3.2485447959084092e+02
         xcrref[5] = 3.2707541254659363e+03

#---------------------------------------------------------------------
#  Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------

         if iotype == 0
           xceref[1] = 5.2969847140936856e+01
           xceref[2] = 4.4632896115670668e+00
           xceref[3] = 1.3122573342210174e+01
           xceref[4] = 1.2006925323559144e+01
           xceref[5] = 1.2459576151035986e+02
         else
           xceref[1] = 0.1477545106464e+03
           xceref[2] = 0.1108895555053e+02
           xceref[3] = 0.3698065590331e+02
           xceref[4] = 0.3310505581440e+02
           xceref[5] = 0.3157928282563e+03
         end

#---------------------------------------------------------------------
#    reference data for 162X162X162 grids after 200 time steps,
#    with DT = 1.0d-04
#---------------------------------------------------------------------
        elseif class == "C"

           dtref = 1.0e-4

#---------------------------------------------------------------------
#  Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
         xcrref[1] = 0.62398116551764615e+04
         xcrref[2] = 0.50793239190423964e+03
         xcrref[3] = 0.15423530093013596e+04
         xcrref[4] = 0.13302387929291190e+04
         xcrref[5] = 0.11604087428436455e+05

#---------------------------------------------------------------------
#  Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------

         if iotype == 0
           xceref[1] = 0.16462008369091265e+03
           xceref[2] = 0.11497107903824313e+02
           xceref[3] = 0.41207446207461508e+02
           xceref[4] = 0.37087651059694167e+02
           xceref[5] = 0.36211053051841265e+03
         else
           xceref[1] = 0.2597156483475e+03
           xceref[2] = 0.1985384289495e+02
           xceref[3] = 0.6517950485788e+02
           xceref[4] = 0.5757235541520e+02
           xceref[5] = 0.5215668188726e+03
         end


#---------------------------------------------------------------------
#    reference data for 408x408x408 grids after 250 time steps,
#    with DT = 0.2d-04
#---------------------------------------------------------------------
        elseif class == "D"

           dtref = 0.2e-4

#---------------------------------------------------------------------
#  Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
         xcrref[1] = 0.2533188551738e+05
         xcrref[2] = 0.2346393716980e+04
         xcrref[3] = 0.6294554366904e+04
         xcrref[4] = 0.5352565376030e+04
         xcrref[5] = 0.3905864038618e+05

#---------------------------------------------------------------------
#  Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------

         if iotype == 0
           xceref[1] = 0.3100009377557e+03
           xceref[2] = 0.2424086324913e+02
           xceref[3] = 0.7782212022645e+02
           xceref[4] = 0.6835623860116e+02
           xceref[5] = 0.6065737200368e+03
         else
           xceref[1] = 0.3813781566713e+03
           xceref[2] = 0.3160872966198e+02
           xceref[3] = 0.9593576357290e+02
           xceref[4] = 0.8363391989815e+02
           xceref[5] = 0.7063466087423e+03
         end


#---------------------------------------------------------------------
#    reference data for 1020x1020x1020 grids after 250 time steps,
#    with DT = 0.4d-05
#---------------------------------------------------------------------
        elseif class == "E"

           dtref = 0.4e-5

#---------------------------------------------------------------------
#  Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
         xcrref[1] = 0.9795372484517e+05
         xcrref[2] = 0.9739814511521e+04
         xcrref[3] = 0.2467606342965e+05
         xcrref[4] = 0.2092419572860e+05
         xcrref[5] = 0.1392138856939e+06

#---------------------------------------------------------------------
#  Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------

         if iotype == 0
           xceref[1] = 0.4327562208414e+03
           xceref[2] = 0.3699051964887e+02
           xceref[3] = 0.1089845040954e+03
           xceref[4] = 0.9462517622043e+02
           xceref[5] = 0.7765512765309e+03
         else
#  wr_interval = 5
           xceref[1] = 0.4729898413058e+03
           xceref[2] = 0.4145899331704e+02
           xceref[3] = 0.1192850917138e+03
           xceref[4] = 0.1032746026932e+03
           xceref[5] = 0.8270322177634e+03
#  wr_interval = 10
#          xceref(1) = 0.4718135916251d+03
#          xceref(2) = 0.4132620259096d+02
#          xceref(3) = 0.1189831133503d+03
#          xceref(4) = 0.1030212798803d+03
#          xceref(5) = 0.8255924078458d+03
         end

#---------------------------------------------------------------------
#    reference data for 2560x2560x2560 grids after 250 time steps,
#    with DT = 0.6d-06
#---------------------------------------------------------------------
        elseif class == "F"

           dtref = 0.6e-6

#---------------------------------------------------------------------
#  Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
         xcrref[1] = 0.4240735175585e+06
         xcrref[2] = 0.4348701133212e+05
         xcrref[3] = 0.1078114688845e+06
         xcrref[4] = 0.9142160938556e+05
         xcrref[5] = 0.5879842143431e+06

#---------------------------------------------------------------------
#  Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------

         if iotype == 0
           xceref[1] = 0.5095577042351e+03
           xceref[2] = 0.4557065541652e+02
           xceref[3] = 0.1286632140581e+03
           xceref[4] = 0.1111419378722e+03
           xceref[5] = 0.8720011709356e+03
         end

        else

           verified = false

        end

#---------------------------------------------------------------------
#    verification test for residuals if gridsize is one of 
#    the defined grid sizes above (class .ne. 'U')
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#    Compute the difference of solution values and the known reference 
#    values.
#---------------------------------------------------------------------
        for m = 1:5

           xcrdif[m] = abs((xcr[m]-xcrref[m])/xcrref[m])
           xcedif[m] = abs((xce[m]-xceref[m])/xceref[m])

        end

#---------------------------------------------------------------------
#    Output the comparison of computed results to known cases.
#---------------------------------------------------------------------

        if class != "U"
           @printf(stdout, " Verification being performed for class %s\n", class)
# 1990      format(' Verification being performed for class ', a)
           @printf(stdout, " accuracy setting for epsilon = %20.13E\n", epsilon)
# 2000      format(' accuracy setting for epsilon = ', E20.13)
           verified = (abs(dt-dtref) <= epsilon)
           if !verified
              class = "U"
              @printf(stdout, " DT does not match the reference value of %15.8E\n", dtref)
# 1000         format(' DT does not match the reference value of ',  		                       E15.8)
           end
        else
           @printf(stdout, " Unknown class\n", )
# 1995      format(' Unknown class')
        end


        if class != "U"
           @printf(stdout, " Comparison of RMS-norms of residual\n", )
        else
           @printf(stdout, " RMS-norms of residual\n", )
        end

# 2001   format(' Comparison of RMS-norms of residual')
# 2005   format(' RMS-norms of residual')
        for m = 1:5
           if class == "U"
              @printf(stdout, "          %2i%20.13E\n", m, xcr[m])
           elseif xcrdif[m] != NaN && xcrdif[m] <= epsilon
              @printf(stdout, "          %2i%20.13E%20.13E%20.13E\n", m, xcr[m], xcrref[m], xcrdif[m])
           else
              verified = false
              @printf(stdout, " FAILURE: %2i%20.13E%20.13E%20.13E\n", m, xcr[m], xcrref[m], xcrdif[m])
           end
        end

        if class != "U"
           @printf(stdout, " Comparison of RMS-norms of solution error\n", )
        else
           @printf(stdout, " RMS-norms of solution error\n", )
        end
# 2002   format(' Comparison of RMS-norms of solution error')
# 2006   format(' RMS-norms of solution error')

        for m = 1:5
           if class == "U"
              @printf(stdout, "          %2i%20.13E\n", m, xce[m])
           elseif xcedif[m] != NaN && xcedif[m] <= epsilon
              @printf(stdout, "          %2i%20.13E%20.13E%20.13E\n", m, xce[m], xceref[m], xcedif[m])
           else
              verified = false
              @printf(stdout, " FAILURE: %2i%20.13E%20.13E%20.13E\n", m, xce[m], xceref[m], xcedif[m])
           end
        end

# 2010   format(' FAILURE: ', i2, E20.13, E20.13, E20.13)
# 2011   format('          ', i2, E20.13, E20.13, E20.13)
# 2015   format('          ', i2, E20.13)

        if class == "U"
           @printf(stdout, " No reference values provided\n", )
           @printf(stdout, " No verification performed\n", )
# 2022      format(' No reference values provided')
# 2023      format(' No verification performed')
        elseif verified
           @printf(stdout, " Verification Successful\n", )
# 2020      format(' Verification Successful')
        else
           @printf(stdout, " Verification failed\n", )
# 2021      format(' Verification failed')
        end

        return nothing


        end
