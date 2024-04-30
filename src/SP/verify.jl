function set_class(no_time_steps, x_zones, y_zones, gx_size, gy_size, gz_size)

        if x_zones == sp_class[CLASS_S].x_zones &&             
           y_zones == sp_class[CLASS_S].y_zones &&
           gx_size == sp_class[CLASS_S].gx_size &&
           gy_size == sp_class[CLASS_S].gy_size &&
           gz_size == sp_class[CLASS_S].gz_size &&
           no_time_steps == sp_class[CLASS_S].niter
             class = CLASS_S
        elseif x_zones == sp_class[CLASS_W].x_zones &&             
               y_zones == sp_class[CLASS_W].y_zones &&
               gx_size == sp_class[CLASS_W].gx_size &&
               gy_size == sp_class[CLASS_W].gy_size &&
               gz_size == sp_class[CLASS_W].gz_size &&
               no_time_steps == sp_class[CLASS_W].niter
             class = CLASS_W
        elseif x_zones == sp_class[CLASS_A].x_zones &&             
               y_zones == sp_class[CLASS_A].y_zones &&
               gx_size == sp_class[CLASS_A].gx_size &&
               gy_size == sp_class[CLASS_A].gy_size &&
               gz_size == sp_class[CLASS_A].gz_size &&
               no_time_steps == sp_class[CLASS_A].niter
             class = CLASS_A
        elseif x_zones == sp_class[CLASS_B].x_zones &&             
               y_zones == sp_class[CLASS_B].y_zones &&
               gx_size == sp_class[CLASS_B].gx_size &&
               gy_size == sp_class[CLASS_B].gy_size &&
               gz_size == sp_class[CLASS_B].gz_size &&
               no_time_steps == sp_class[CLASS_B].niter
             class = CLASS_B
        elseif x_zones == sp_class[CLASS_C].x_zones &&             
               y_zones == sp_class[CLASS_C].y_zones &&
               gx_size == sp_class[CLASS_C].gx_size &&
               gy_size == sp_class[CLASS_C].gy_size &&
               gz_size == sp_class[CLASS_C].gz_size &&
               no_time_steps == sp_class[CLASS_C].niter
             class = CLASS_C
        elseif x_zones == sp_class[CLASS_D].x_zones &&             
               y_zones == sp_class[CLASS_D].y_zones &&
               gx_size == sp_class[CLASS_D].gx_size &&
               gy_size == sp_class[CLASS_D].gy_size &&
               gz_size == sp_class[CLASS_D].gz_size &&
               no_time_steps == sp_class[CLASS_D].niter
             class = CLASS_D
        elseif x_zones == sp_class[CLASS_E].x_zones &&             
               y_zones == sp_class[CLASS_E].y_zones &&
               gx_size == sp_class[CLASS_E].gx_size &&
               gy_size == sp_class[CLASS_E].gy_size &&
               gz_size == sp_class[CLASS_E].gz_size &&
               no_time_steps == sp_class[CLASS_E].niter
             class = CLASS_E
        elseif x_zones == sp_class[CLASS_F].x_zones &&             
               y_zones == sp_class[CLASS_F].y_zones &&
               gx_size == sp_class[CLASS_F].gx_size &&
               gy_size == sp_class[CLASS_F].gy_size &&
               gz_size == sp_class[CLASS_F].gz_size &&
               no_time_steps == sp_class[CLASS_F].niter
             class = CLASS_F
        else
             class = CLASS_UNDEFINED
        end

        return class

end

# AT CLUSTER MASTER
function putVerificationReport(xcr_node, xce_node)
   remotecall(putVerificationReport, 1, clusterid, xcr_node, xce_node; role=:worker)
end

# AT DRIVER MASTER
function putVerificationReport(clusterid, xcr_cluster, xce_cluster)
   lock(verified_signal[clusterid+1])
   try
      lock(lk_update_verify)
      try
         xcr .= xcr + xcr_cluster
         xce .= xce + xce_cluster
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
#---------------------------------------------------------------------
# AT DRIVER MASTER
function verify(class, dt, no_time_steps)

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

        for m = 1:5
           xcrref[m] = 1.0
           xceref[m] = 1.0
        end

#---------------------------------------------------------------------
#    reference data for class S
#---------------------------------------------------------------------
        if class == CLASS_S
           dtref = 1.5e-2
           niterref = 100

#---------------------------------------------------------------------
#    Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.7698876173566e+01
           xcrref[2] = 0.1517766790280e+01
           xcrref[3] = 0.2686805141546e+01
           xcrref[4] = 0.1893688083690e+01
           xcrref[5] = 0.1369739859738e+02

#---------------------------------------------------------------------
#    Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.9566808043467e+01
           xceref[2] = 0.3894109553741e+01
           xceref[3] = 0.4516022447464e+01
           xceref[4] = 0.4099103995615e+01
           xceref[5] = 0.7776038881521e+01

#---------------------------------------------------------------------
#    reference data for class W
#---------------------------------------------------------------------
        elseif class == CLASS_W
           dtref = 1.5e-3
           niterref = 400

#---------------------------------------------------------------------
#    Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.1887636218359e+03
           xcrref[2] = 0.1489637963542e+02
           xcrref[3] = 0.4851711701400e+02
           xcrref[4] = 0.3384633608154e+02
           xcrref[5] = 0.4036632495857e+03

#---------------------------------------------------------------------
#    Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.2975895149929e+02
           xceref[2] = 0.1341508175806e+02
           xceref[3] = 0.1585310846491e+02
           xceref[4] = 0.1450916426713e+02
           xceref[5] = 0.5854137431023e+02

#---------------------------------------------------------------------
#    reference data for class A
#---------------------------------------------------------------------
        elseif class == CLASS_A
           dtref = 1.5e-3
           niterref = 400

#---------------------------------------------------------------------
#    Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.2800097900548e+03
           xcrref[2] = 0.2268349014438e+02
           xcrref[3] = 0.7000852739901e+02
           xcrref[4] = 0.5000771004061e+02
           xcrref[5] = 0.5552068537578e+03

#---------------------------------------------------------------------
#    Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.3112046666578e+02
           xceref[2] = 0.1172197785348e+02
           xceref[3] = 0.1486616708032e+02
           xceref[4] = 0.1313680576292e+02
           xceref[5] = 0.7365834058154e+02

#---------------------------------------------------------------------
#    reference data for class B
#---------------------------------------------------------------------
        elseif class == CLASS_B
           dtref = 1.0e-3
           niterref = 400

#---------------------------------------------------------------------
#    Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.5190422977921e+04
           xcrref[2] = 0.3655458539065e+03
           xcrref[3] = 0.1261126592633e+04
           xcrref[4] = 0.1002038338842e+04
           xcrref[5] = 0.1075902511165e+05

#---------------------------------------------------------------------
#    Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.5469182054223e+03
           xceref[2] = 0.4983658028989e+02
           xceref[3] = 0.1418301776602e+03
           xceref[4] = 0.1097717156175e+03
           xceref[5] = 0.1260195162174e+04

#---------------------------------------------------------------------
#    reference data for class C
#---------------------------------------------------------------------
        elseif class == CLASS_C
           dtref = 0.67e-3
           niterref = 400

#---------------------------------------------------------------------
#    Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.5886814493676e+05
           xcrref[2] = 0.3967324375474e+04
           xcrref[3] = 0.1444126529019e+05
           xcrref[4] = 0.1210582211196e+05
           xcrref[5] = 0.1278941567976e+06

#---------------------------------------------------------------------
#    Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.6414069213021e+04
           xceref[2] = 0.4069468353404e+03
           xceref[3] = 0.1585311908719e+04
           xceref[4] = 0.1270243185759e+04
           xceref[5] = 0.1441398372869e+05

#---------------------------------------------------------------------
#    reference data for class D
#---------------------------------------------------------------------
        elseif class == CLASS_D
           dtref = 0.3e-3
           niterref = 500

#---------------------------------------------------------------------
#    Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.7650595424723e+06
           xcrref[2] = 0.5111519817683e+05
           xcrref[3] = 0.1857213937602e+06
           xcrref[4] = 0.1624096784059e+06
           xcrref[5] = 0.1642416844328e+07

#---------------------------------------------------------------------
#    Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.8169589578340e+05
           xceref[2] = 0.5252150843148e+04
           xceref[3] = 0.1984739188642e+05
           xceref[4] = 0.1662852404547e+05
           xceref[5] = 0.1761381855235e+06

#---------------------------------------------------------------------
#    reference data for class E
#---------------------------------------------------------------------
        elseif class == CLASS_E
           dtref = 0.2e-3
           niterref = 500

#---------------------------------------------------------------------
#    Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.5058298119039e+07
           xcrref[2] = 0.3576837494299e+06
           xcrref[3] = 0.1230856227329e+07
           xcrref[4] = 0.1093895671677e+07
           xcrref[5] = 0.1073671658903e+08

#---------------------------------------------------------------------
#    Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.5288293042051e+06
           xceref[2] = 0.3471875724140e+05
           xceref[3] = 0.1282998930808e+06
           xceref[4] = 0.1095483394612e+06
           xceref[5] = 0.1129716454231e+07

#---------------------------------------------------------------------
#    reference data for class F
#---------------------------------------------------------------------
        elseif class == CLASS_F
           dtref = 0.1e-3
           niterref = 500

#---------------------------------------------------------------------
#    Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.3974469160412e+08
           xcrref[2] = 0.3260760921834e+07
           xcrref[3] = 0.9756215393494e+07
           xcrref[4] = 0.8278472138497e+07
           xcrref[5] = 0.7547269314441e+08

#---------------------------------------------------------------------
#    Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.3475757666334e+07
           xceref[2] = 0.2386799228183e+06
           xceref[3] = 0.8436705443034e+06
           xceref[4] = 0.7339112115118e+06
           xceref[5] = 0.7327832757877e+07

           if no_time_steps == 50

               niterref = 50
               xcrref[1] = 0.3198801286787e+09
               xcrref[2] = 0.3435698123358e+08
               xcrref[3] = 0.8489831174901e+08
               xcrref[4] = 0.6940707552477e+08
               xcrref[5] = 0.4478684103255e+09

               xceref[1] = 0.6761099692230e+07
               xceref[2] = 0.5361561494769e+06
               xceref[3] = 0.1662878706114e+07
               xceref[4] = 0.1443852092060e+07
               xceref[5] = 0.1260678700480e+08

           end
        else
           dtref = 0.0e0
           niterref = 0
           verified = false
        end

#---------------------------------------------------------------------
#    Compute the difference of solution values and the known reference values.
#---------------------------------------------------------------------

        # VETORIZED
        xcrdif = abs.((xcr .- xcrref)./xcrref)
        xcedif = abs.((xce .- xceref)./xceref)

#---------------------------------------------------------------------
#    Output the comparison of computed results to known cases.
#---------------------------------------------------------------------

        @printf(stdout, " Verification being performed for class %s\n", class)
        @printf(stdout, " accuracy setting for epsilon = %20.13E\n", epsilon)
        if abs(dt-dtref) > epsilon
           verified = false
           @printf(stdout, " DT does not match the reference value of %15.8E\n", dtref)
        elseif no_time_steps != niterref
           verified = false
           @printf(stdout, " NITER does not match the reference value of %5i\n", niterref)
        end

        @printf(stdout, " Comparison of RMS-norms of residual\n", )

        for m = 1:5
           if xcrdif[m] != NaN && xcrdif[m] <= epsilon
              @printf(stdout, "          %2i%20.13E%20.13E%20.13E\n", m, xcr[m], xcrref[m], xcrdif[m])
           else
              verified = false
              @printf(stdout, " FAILURE: %2i%20.13E%20.13E%20.13E\n", m, xcr[m], xcrref[m], xcrdif[m])
           end
        end

        @printf(stdout, " Comparison of RMS-norms of solution error\n", )

        for m = 1:5
           if  xcedif[m] != NaN && xcedif[m] <= epsilon
              @printf(stdout, "          %2i%20.13E%20.13E%20.13E\n", m, xce[m], xceref[m], xcedif[m])
           else
              verified = false
              @printf(stdout, " FAILURE: %2i%20.13E%20.13E%20.13E\n", m, xce[m], xceref[m], xcedif[m])
           end
        end

        if verified
           @printf(stdout, " Verification Successful\n", )
        else
           @printf(stdout, " Verification failed\n", )
        end

        return verified
end


# AT NODE WORKER
function verify(dt, ss, sr, b_size, 
                proc_num_zones, proc_zone_id,
                rho_i, us, vs, ws, speed, qs, square,
                rhs, forcing, u, nx, ny, nz)

        xce = Array{FloatType}(undef, 5)
        xcr = Array{FloatType}(undef, 5)
        xce_sub = Array{FloatType}(undef, 5)
        xcr_sub = Array{FloatType}(undef, 5)

#---------------------------------------------------------------------
#   compute the error norm and the residual norm
#---------------------------------------------------------------------

        for m = 1:5
          xcr[m] = 0.0e0
          xce[m] = 0.0e0
        end

        for iz = 1:proc_num_zones          

          zone = proc_zone_id[iz]

          error_norm(iz, xce_sub, [nx[zone], ny[zone], nz[zone]])
          
          copy_faces(false, iz,
                  Val(no_nodes),
                  Val(ncells),
                  successor[iz],
                  predecessor[iz],
                  cell_size[iz],
                  cell_start[iz],
                  cell_end[iz],
                  cell_coord[iz],
                  cell_low[iz],
                  cell_high[iz],
                  u[iz],
                  rhs[iz],
                  rho_i[iz],
                  us[iz],
                  vs[iz],
                  ws[iz],
                  square[iz],
                  qs[iz],
                  ainv[iz],
                  speed[iz],
                  forcing[iz],
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
                  in_buffer[iz],
                  out_buffer[iz],
                  Array{MPI.Request}(undef,12),
                  timeron,
                  comm_rhs[iz],
                  ss[iz],
                  sr[iz],
                  b_size[iz],
               )

          rhs_norm(iz, xcr_sub, [nx[zone], ny[zone], nz[zone]])

          for m = 1:5
            xcr[m] = xcr[m] + xcr_sub[m] / dt
            xce[m] = xce[m] + xce_sub[m]
          end

        end

        if node == root
           remotecall(putVerificationReport, 1, xcr, xce; role = :worker)
        end
end
