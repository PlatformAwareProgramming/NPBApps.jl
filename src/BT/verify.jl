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
           
   xcrref = Array{Float64}(undef, 5)
   xceref = Array{Float64}(undef, 5)
   xcrdif = Array{Float64}(undef, 5)
   xcedif = Array{Float64}(undef, 5)

#---------------------------------------------------------------------
#   tolerance level
#---------------------------------------------------------------------
   epsilon = 1.0e-08

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
      dtref = 0.010e0
      niterref = 60

#---------------------------------------------------------------------
#  Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.1047687395830e+04
           xcrref[2] = 0.9419911314792e+02
           xcrref[3] = 0.2124737403068e+03
           xcrref[4] = 0.1422173591794e+03
           xcrref[5] = 0.1135441572375e+04

#---------------------------------------------------------------------
#  Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.1775416062982e+03
           xceref[2] = 0.1875540250835e+02
           xceref[3] = 0.3863334844506e+02
           xceref[4] = 0.2634713890362e+02
           xceref[5] = 0.1965566269675e+03


#---------------------------------------------------------------------
#    reference data for class W
#---------------------------------------------------------------------
   elseif class == CLASS_W
      dtref = 0.0008e0
      niterref = 200

#---------------------------------------------------------------------
#  Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.5562611195402e+05
           xcrref[2] = 0.5151404119932e+04
           xcrref[3] = 0.1080453907954e+05
           xcrref[4] = 0.6576058591929e+04
           xcrref[5] = 0.4528609293561e+05

#---------------------------------------------------------------------
#  Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.7185154786403e+04
           xceref[2] = 0.7040472738068e+03
           xceref[3] = 0.1437035074443e+04
           xceref[4] = 0.8570666307849e+03
           xceref[5] = 0.5991235147368e+04


#---------------------------------------------------------------------
#    reference data for class A
#---------------------------------------------------------------------
   elseif class == CLASS_A
      dtref = 0.0008e0
      niterref = 200

#---------------------------------------------------------------------
#  Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.5536703889522e+05
           xcrref[2] = 0.5077835038405e+04
           xcrref[3] = 0.1067391361067e+05
           xcrref[4] = 0.6441179694972e+04
           xcrref[5] = 0.4371926324069e+05

#---------------------------------------------------------------------
#  Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.6716797714343e+04
           xceref[2] = 0.6512687902160e+03
           xceref[3] = 0.1332930740128e+04
           xceref[4] = 0.7848302089180e+03
           xceref[5] = 0.5429053878818e+04

#---------------------------------------------------------------------
#    reference data for class B
#---------------------------------------------------------------------
   elseif class == CLASS_B
      dtref = 0.0003e0
      niterref = 200

#---------------------------------------------------------------------
#  Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.4461388343844e+06
           xcrref[2] = 0.3799759138035e+05
           xcrref[3] = 0.8383296623970e+05
           xcrref[4] = 0.5301970201273e+05
           xcrref[5] = 0.3618106851311e+06

#---------------------------------------------------------------------
#  Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.4496733567600e+05
           xceref[2] = 0.3892068540524e+04
           xceref[3] = 0.8763825844217e+04
           xceref[4] = 0.5599040091792e+04
           xceref[5] = 0.4082652045598e+05

#---------------------------------------------------------------------
#    reference data for class C
#---------------------------------------------------------------------
   elseif class == CLASS_C
      dtref = 0.0001e0
      niterref = 200

#---------------------------------------------------------------------
#  Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.3457703287806e+07
           xcrref[2] = 0.3213621375929e+06
           xcrref[3] = 0.7002579656870e+06
           xcrref[4] = 0.4517459627471e+06
           xcrref[5] = 0.2818715870791e+07

#---------------------------------------------------------------------
#  Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.2059106993570e+06
           xceref[2] = 0.1680761129461e+05
           xceref[3] = 0.4080731640795e+05
           xceref[4] = 0.2836541076778e+05
           xceref[5] = 0.2136807610771e+06

#---------------------------------------------------------------------
#    reference data for class D
#---------------------------------------------------------------------
   elseif class == "D"
      dtref =  0.00002e0
      niterref = 250

#---------------------------------------------------------------------
#  Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.4250417034981e+08
           xcrref[2] = 0.4293882192175e+07
           xcrref[3] = 0.9121841878270e+07
           xcrref[4] = 0.6201357771439e+07
           xcrref[5] = 0.3474801891304e+08

#---------------------------------------------------------------------
#  Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.9462418484583e+06
           xceref[2] = 0.7884728947105e+05
           xceref[3] = 0.1902874461259e+06
           xceref[4] = 0.1361858029909e+06
           xceref[5] = 0.9816489456253e+06


#---------------------------------------------------------------------
#    reference data for class E
#---------------------------------------------------------------------
   elseif class == CLASS_E
      dtref = 0.000004e0
      niterref = 250

#---------------------------------------------------------------------
#  Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.5744815962469e+09
           xcrref[2] = 0.6088696479719e+08
           xcrref[3] = 0.1276325224438e+09
           xcrref[4] = 0.8947040105616e+08
           xcrref[5] = 0.4726115284807e+09

#---------------------------------------------------------------------
#  Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.4114447054461e+07
           xceref[2] = 0.3570776728190e+06
           xceref[3] = 0.8465106191458e+06
           xceref[4] = 0.6147182273817e+06
           xceref[5] = 0.4238908025163e+07

#---------------------------------------------------------------------
#    reference data for class F
#---------------------------------------------------------------------
   elseif class == CLASS_F
      dtref = 0.000001e0
      niterref = 250

#---------------------------------------------------------------------
#  Reference values of RMS-norms of residual.
#---------------------------------------------------------------------
           xcrref[1] = 0.6524078317845e+10
           xcrref[2] = 0.7020439279514e+09
           xcrref[3] = 0.1467588422194e+10
           xcrref[4] = 0.1042973064137e+10
           xcrref[5] = 0.5411102201141e+10

#---------------------------------------------------------------------
#  Reference values of RMS-norms of solution error.
#---------------------------------------------------------------------
           xceref[1] = 0.1708795375347e+08
           xceref[2] = 0.1514359936802e+07
           xceref[3] = 0.3552878359250e+07
           xceref[4] = 0.2594549582184e+07
           xceref[5] = 0.1749809607845e+08

      if no_time_steps == 50
           niterref = 250
           xcrref[1] = 0.3565049484400e+11
           xcrref[2] = 0.3752029586145e+10
           xcrref[3] = 0.7805935552197e+10
           xcrref[4] = 0.5685995438056e+10
           xcrref[5] = 0.2908811276266e+11

           xceref[1] = 0.1805995755490e+08
           xceref[2] = 0.1632306899424e+07
           xceref[3] = 0.3778610439036e+07
           xceref[4] = 0.2749319818549e+07
           xceref[5] = 0.1814401049296e+08
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
           rho_i, us, vs, ws, qs, square,
           rhs, forcing, u, nx, ny, nz, requests)

   xce = Array{Float64}(undef, 5)
   xcr = Array{Float64}(undef, 5)
   xce_sub = Array{Float64}(undef, 5)
   xcr_sub = Array{Float64}(undef, 5)

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
     
       copy_faces(ss[iz], 
                  sr[iz], 
                  b_size[iz],
                  cell_coord[iz],
                  cell_size[iz],
                  cell_start[iz],
                  cell_end[iz],
                  forcing[iz],        
                  u[iz],
                  rhs[iz],
                  in_buffer[iz],
                  out_buffer[iz],
                  us[iz],
                  vs[iz],
                  ws[iz],
                  qs[iz],
                  rho_i[iz],
                  square[iz],
                  timeron,
                  dt,
                  Val(ncells),
                  tx2,
                  ty2,
                  tz2,
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
                  Val(no_nodes), 
                  comm_rhs[iz],
                  predecessor[iz],
                  successor[iz],
                  requests[iz],
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





























