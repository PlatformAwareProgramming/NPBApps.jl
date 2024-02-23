function zone_setup(x_zones, y_zones, gx_size, gy_size, gz_size, nx, nxmax, ny, nz, ratio, npb_verbose)

#=       if abs(ratio-1.0e0) > 1.0e-10

#        compute zone stretching only if the prescribed zone size ratio is substantially larger than unity       

         x_r   = exp(log(ratio)/(x_zones-1))
         y_r   = exp(log(ratio)/(y_zones-1))
         x_smallest = float(gx_size)*(x_r-1.0e0)/(x_r^x_zones-1.0e0)
         y_smallest = float(gy_size)*(y_r-1.0e0)/(y_r^y_zones-1.0e0)

#        compute tops of intervals, using a slightly tricked rounding to make sure that the intervals are increasing monotonically in size

         for i = 1:x_zones
            x_end[i] = round(x_smallest*(x_r^i-1.0e0)/(x_r-1.0e0)+0.45e0)
         end

         for j = 1:y_zones
            y_end[j] = round(y_smallest*(y_r^j-1.0e0)/(y_r-1.0e0)+0.45e0) 
         end

         @info "x_end=$x_end, y_end=$y_end"

       else=#

#        compute essentially equal sized zone dimensions

         for i = 1:x_zones
           x_end[i]   = (i*gx_size)/x_zones
         end

         for j = 1:y_zones
           y_end[j]   = (j*gy_size)/y_zones
         end

       #end

       x_start[1] = 1
       for i = 1:x_zones
          if (i != x_zones) x_start[i+1] = x_end[i] + 1 end
          x_size[i]  = x_end[i] - x_start[i] + 1
       end

       y_start[1] = 1
       for j = 1:y_zones
          if (j != y_zones) y_start[j+1] = y_end[j] + 1 end
          y_size[j] = y_end[j] - y_start[j] + 1
       end

       if npb_verbose > 1 
          @printf(stdout, "\n Zone sizes:\n", )
       end

       for j = 1:y_zones
         for i = 1:x_zones
           zone_no = (i-1)+(j-1)*x_zones+1
           nx[zone_no] = x_size[i]
           nxmax[zone_no] = nx[zone_no] + 1 - mod(nx[zone_no], 2)
           ny[zone_no] = y_size[j]
           nz[zone_no] = gz_size

           id_west  = mod(i-2+x_zones, x_zones)
           id_east  = mod(i,          x_zones)
           jd_south = mod(j-2+y_zones, y_zones)
           jd_north = mod(j,          y_zones)
           iz_west[zone_no] = id_west +  (j-1)*x_zones + 1
           iz_east[zone_no] = id_east +  (j-1)*x_zones + 1
           iz_south[zone_no] = (i-1) + jd_south*x_zones + 1
           iz_north[zone_no] = (i-1) + jd_north*x_zones + 1

           if npb_verbose > 1
             @printf(stdout, "%5i:  %5i x%5i x%5i\n", zone_no, nx[zone_no], ny[zone_no], nz[zone_no])
           end
         end
       end

       return nothing
end



