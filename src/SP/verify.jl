#---------------------------------------------------------------------
#  set problem class based on problem size
#---------------------------------------------------------------------

function set_class(no_time_steps, grid_points)

        if (grid_points[1]  == 12     ) &&(
             grid_points[2]  == 12     ) &&(
             grid_points[3]  == 12     ) &&(
             no_time_steps   == 100    )

           class = CLASS_S

        elseif (grid_points[1] == 36) &&(
                 grid_points[2] == 36) &&(
                 grid_points[3] == 36) &&(
                 no_time_steps  == 400)

           class = CLASS_W

        elseif (grid_points[1] == 64) &&(
                 grid_points[2] == 64) &&(
                 grid_points[3] == 64) &&(
                 no_time_steps  == 400)

           class = CLASS_A

        elseif (grid_points[1] == 102) &&(
                 grid_points[2] == 102) &&(
                 grid_points[3] == 102) &&(
                 no_time_steps  == 400)

           class = CLASS_B

        elseif (grid_points[1] == 162) &&(
                 grid_points[2] == 162) &&(
                 grid_points[3] == 162) &&(
                 no_time_steps  == 400)

           class = CLASS_C

        elseif (grid_points[1] == 408) &&(
                 grid_points[2] == 408) &&(
                 grid_points[3] == 408) &&(
                 no_time_steps  == 500)

           class = CLASS_D

        elseif (grid_points[1] == 1020) &&(
                 grid_points[2] == 1020) &&(
                 grid_points[3] == 1020) &&(
                 no_time_steps  == 500)

           class = CLASS_E

        elseif (grid_points[1] == 2560) &&(
                 grid_points[2] == 2560) &&(
                 grid_points[3] == 2560) &&(
                 no_time_steps  == 500)

           class = CLASS_F

        else

           class = CLASS_UNDEFINED

        end

        return class

end

#---------------------------------------------------------------------
#  verification routine                         
#---------------------------------------------------------------------

function verify(no_nodes, node,
                ncells, 
                class, 
                grid_points,
                cell_coord, cell_low, cell_high, cell_start, cell_end, cell_size,
                u,
                rhs,
                rho_i,
                us,
                vs,
                ws,
                square,
                qs,
                ainv,
                speed,
                forcing,
                in_buffer,
                out_buffer,
                dt,
                ss,
                sr,
                b_size,
                tx2, ty2, tz2,
                c1, c2, c1c2,
                dx1tx1, dx2tx1, dx3tx1, dx4tx1, dx5tx1,
                dy1ty1, dy2ty1, dy3ty1, dy4ty1, dy5ty1,
                dz1tz1, dz2tz1, dz3tz1, dz4tz1, dz5tz1,
                dnxm1, dnym1, dnzm1,
                xxcon2, xxcon3, xxcon4, xxcon5,
                yycon2, yycon3, yycon4, yycon5,
                zzcon2, zzcon3, zzcon4, zzcon5,
                dssp,
                con43,
                timeron,
                comm_setup,
                comm_rhs
)

        xcrref = Array{Float64}(undef, 5)
        xceref = Array{Float64}(undef, 5)
        xcrdif = Array{Float64}(undef, 5)
        xcedif = Array{Float64}(undef, 5)
        xce = Array{Float64}(undef, 5)
        xcr = Array{Float64}(undef, 5)

#---------------------------------------------------------------------
#   tolerance level
#---------------------------------------------------------------------
        epsilon = 1.0e-08

#---------------------------------------------------------------------
#   compute the error norm and the residual norm, and exit if not printing
#---------------------------------------------------------------------
        error_norm(ncells, u, xce, grid_points, cell_low, cell_high, dnxm1, dnym1, dnzm1, comm_setup)
        
        copy_faces(Val(no_nodes),
                  Val(ncells),
                  successor, # ::Vector{Int64},
                  predecessor, # ::Vector{Int64},
                  cell_size,
                  cell_start,
                  cell_end,
                  cell_coord,
                  u,
                  rhs,
                  rho_i,
                  us,
                  vs,
                  ws,
                  square,
                  qs,
                  ainv,
                  speed,
                  forcing,
                  dt,
                  tx2,
                  ty2,
                  tz2,
                  c1,
                  c2,
                  c1c2,
                  dx1tx1,
                  dx2tx1,
                  dx3tx1, 
                  dx4tx1,
                  dx5tx1,
                  dy1ty1,
                  dy2ty1,
                  dy3ty1,
                  dy4ty1,
                  dy5ty1,
                  dz1tz1,
                  dz2tz1,
                  dz3tz1,
                  dz4tz1,
                  dz5tz1,
                  xxcon2,
                  xxcon3,
                  xxcon4,
                  xxcon5,
                  yycon2,
                  yycon3,
                  yycon4,
                  yycon5,
                  zzcon2,
                  zzcon3,
                  zzcon4,
                  zzcon5,
                  dssp,
                  con43,
                  in_buffer,
                  out_buffer,
                  Array{MPI.Request}(undef,12),
                  timeron,
                  comm_rhs,
                  ss,
                  sr,
                  b_size,
               )

        rhs_norm(ncells, rhs, xcr, grid_points, cell_start, cell_end, cell_size, comm_setup)

        # VECTORIZED
        xcr .= xcr ./ dt

        if (node != 0) return end

        verified = true
       
        xcrref .= 1.0
        xceref .= 1.0

#---------------------------------------------------------------------
#    reference data for 12X12X12 grids after 100 time steps, with DT = 1.50d-02
#---------------------------------------------------------------------
        if class == CLASS_S

           dtref = 1.5e-2

#---------------------------------------------------------------------
#    Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 2.7470315451339479e-02
           xcrref[2] = 1.0360746705285417e-02
           xcrref[3] = 1.6235745065095532e-02
           xcrref[4] = 1.5840557224455615e-02
           xcrref[5] = 3.4849040609362460e-02

#---------------------------------------------------------------------
#    Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 2.7289258557377227e-05
           xceref[2] = 1.0364446640837285e-05
           xceref[3] = 1.6154798287166471e-05
           xceref[4] = 1.5750704994480102e-05
           xceref[5] = 3.4177666183390531e-05

#---------------------------------------------------------------------
#    reference data for 36X36X36 grids after 400 time steps, with DT = 1.5d-03
#---------------------------------------------------------------------
        elseif class == CLASS_W

           dtref = 1.5e-3

#---------------------------------------------------------------------
#    Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.1893253733584e-02
           xcrref[2] = 0.1717075447775e-03
           xcrref[3] = 0.2778153350936e-03
           xcrref[4] = 0.2887475409984e-03
           xcrref[5] = 0.3143611161242e-02

#---------------------------------------------------------------------
#    Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.7542088599534e-04
           xceref[2] = 0.6512852253086e-05
           xceref[3] = 0.1049092285688e-04
           xceref[4] = 0.1128838671535e-04
           xceref[5] = 0.1212845639773e-03

#---------------------------------------------------------------------
#    reference data for 64X64X64 grids after 400 time steps, with DT = 1.5d-03
#---------------------------------------------------------------------
        elseif class == CLASS_A

           dtref = 1.5e-3

#---------------------------------------------------------------------
#    Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 2.4799822399300195e0
           xcrref[2] = 1.1276337964368832e0
           xcrref[3] = 1.5028977888770491e0
           xcrref[4] = 1.4217816211695179e0
           xcrref[5] = 2.1292113035138280e0

#---------------------------------------------------------------------
#    Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 1.0900140297820550e-04
           xceref[2] = 3.7343951769282091e-05
           xceref[3] = 5.0092785406541633e-05
           xceref[4] = 4.7671093939528255e-05
           xceref[5] = 1.3621613399213001e-04

#---------------------------------------------------------------------
#    reference data for 102X102X102 grids after 400 time steps,
#    with DT = 1.0d-03
#---------------------------------------------------------------------
        elseif class == CLASS_B

           dtref = 1.0e-3

#---------------------------------------------------------------------
#    Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.6903293579998e+02
           xcrref[2] = 0.3095134488084e+02
           xcrref[3] = 0.4103336647017e+02
           xcrref[4] = 0.3864769009604e+02
           xcrref[5] = 0.5643482272596e+02

#---------------------------------------------------------------------
#    Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.9810006190188e-02
           xceref[2] = 0.1022827905670e-02
           xceref[3] = 0.1720597911692e-02
           xceref[4] = 0.1694479428231e-02
           xceref[5] = 0.1847456263981e-01

#---------------------------------------------------------------------
#    reference data for 162X162X162 grids after 400 time steps,
#    with DT = 0.67d-03
#---------------------------------------------------------------------
        elseif class == CLASS_C

           dtref = 0.67e-3

#---------------------------------------------------------------------
#    Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.5881691581829e+03
           xcrref[2] = 0.2454417603569e+03
           xcrref[3] = 0.3293829191851e+03
           xcrref[4] = 0.3081924971891e+03
           xcrref[5] = 0.4597223799176e+03

#---------------------------------------------------------------------
#    Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.2598120500183e+00
           xceref[2] = 0.2590888922315e-01
           xceref[3] = 0.5132886416320e-01
           xceref[4] = 0.4806073419454e-01
           xceref[5] = 0.5483377491301e+00

#---------------------------------------------------------------------
#    reference data for 408X408X408 grids after 500 time steps,
#    with DT = 0.3d-03
#---------------------------------------------------------------------
        elseif class == CLASS_D

           dtref = 0.30e-3

#---------------------------------------------------------------------
#    Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.1044696216887e+05
           xcrref[2] = 0.3204427762578e+04
           xcrref[3] = 0.4648680733032e+04
           xcrref[4] = 0.4238923283697e+04
           xcrref[5] = 0.7588412036136e+04

#---------------------------------------------------------------------
#    Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.5089471423669e+01
           xceref[2] = 0.5323514855894e+00
           xceref[3] = 0.1187051008971e+01
           xceref[4] = 0.1083734951938e+01
           xceref[5] = 0.1164108338568e+02

#---------------------------------------------------------------------
#    reference data for 1020X1020X1020 grids after 500 time steps,
#    with DT = 0.1d-03
#---------------------------------------------------------------------
        elseif class == CLASS_E

           dtref = 0.10e-3

#---------------------------------------------------------------------
#    Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.6255387422609e+05
           xcrref[2] = 0.1495317020012e+05
           xcrref[3] = 0.2347595750586e+05
           xcrref[4] = 0.2091099783534e+05
           xcrref[5] = 0.4770412841218e+05

#---------------------------------------------------------------------
#    Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.6742735164909e+02
           xceref[2] = 0.5390656036938e+01
           xceref[3] = 0.1680647196477e+02
           xceref[4] = 0.1536963126457e+02
           xceref[5] = 0.1575330146156e+03

#---------------------------------------------------------------------
#    reference data for 2560X2560X2560 grids after 500 time steps,
#    with DT = 0.1d-03
#---------------------------------------------------------------------
        elseif class == CLASS_F

           dtref = 0.15e-4

#---------------------------------------------------------------------
#    Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.9281628449462e+05
           xcrref[2] = 0.2230152287675e+05
           xcrref[3] = 0.3493102358632e+05
           xcrref[4] = 0.3114096186689e+05
           xcrref[5] = 0.7424426448298e+05

#---------------------------------------------------------------------
#    Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.2683717702444e+03
           xceref[2] = 0.2030647554028e+02
           xceref[3] = 0.6734864248234e+02
           xceref[4] = 0.5947451301640e+02
           xceref[5] = 0.5417636652565e+03

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
        #=for m = 1:5

           xcrdif[m] = abs((xcr[m]-xcrref[m])/xcrref[m])
           xcedif[m] = abs((xce[m]-xceref[m])/xceref[m])

        end=#
 
        # VETORIZED
        xcrdif = abs.((xcr .- xcrref)./xcrref)
        xcedif = abs.((xce .- xceref)./xceref)

#---------------------------------------------------------------------
#    Output the comparison of computed results to known cases.
#---------------------------------------------------------------------

        if class != CLASS_UNDEFINED
           @printf(stdout, " Verification being performed for class %s\n", class)
           @printf(stdout, " accuracy setting for epsilon = %20.13E\n", epsilon)
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
              @printf(stdout, "          %2i%20.13E\n", m, xcr[m])
           elseif xcrdif[m] != NaN && xcrdif[m] <= epsilon
              @printf(stdout, "          %2i%20.13E%20.13E%20.13E\n", m, xcr[m], xcrref[m], xcrdif[m])
           else
              verified = false
              @printf(stdout, " FAILURE: %2i%20.13E%20.13E%20.13E\n", m, xcr[m], xcrref[m], xcrdif[m])
           end
        end

        if class != CLASS_UNDEFINED
           @printf(stdout, " Comparison of RMS-norms of solution error\n", )
        else
           @printf(stdout, " RMS-norms of solution error\n", )
        end

        for m = 1:5
           if class == CLASS_UNDEFINED
              @printf(stdout, "          %2i%20.13E\n", m, xce[m])
           elseif xcedif[m] != NaN && xcedif[m] <= epsilon
              @printf(stdout, "          %2i%20.13E%20.13E%20.13E\n", m, xce[m], xceref[m], xcedif[m])
           else
              verified = false
              @printf(stdout, " FAILURE: %2i%20.13E%20.13E%20.13E\n", m, xce[m], xceref[m], xcedif[m])
           end
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
