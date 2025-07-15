#---------------------------------------------------------------------
#
#   initialize MPI and establish rank and size
#
# This is a module in the MPI implementation of LUSSOR
# pseudo application from the NAS Parallel Benchmarks. 
#
#---------------------------------------------------------------------

 function init_comm(proc_num_zones)

#---------------------------------------------------------------------
#    initialize MPI communication
#---------------------------------------------------------------------
      MPI.Init_thread(MPI.THREAD_MULTIPLE)

#---------------------------------------------------------------------
#     get a process grid that requires a (nx*ny) number of procs.
#     excess ranks are marked as inactive.
#---------------------------------------------------------------------

      global xdim, ydim, no_nodes, total_nodes, node, comm_setup, active = get_active_nprocs(MPI.COMM_WORLD, 2)

      if (!active) return end

      global comm_exch = MPI.Comm_dup(comm_setup)
      global comm_solve = Array{MPI.Comm}(undef, proc_num_zones)
      for iz = 1:proc_num_zones
          comm_solve[iz] = MPI.Comm_dup(comm_setup)
      end

#---------------------------------------------------------------------
#   establish the global rank of this process and the group size
#---------------------------------------------------------------------
      global id = node
      global num = no_nodes
      global root = 0

      global ndim   = nodedim(num)

 #=     if !convertdouble
         dp_type = MPI_DOUBLE_PRECISION
      else
         dp_type = MPI_REAL
      end
=#

      return nothing
end

function get_comm_index(zone, iproc, x_zones, y_zones, zone_proc_id, x_size, y_size)
   #
      #      use sp_data
      #      use mpinpb
      #
      #      implicit none
      #
      #  Calculate the communication index of a zone within a processor group
      #
      #      integer zone, iproc, comm_index
      #
      #     local variables
      #      integer izone, jzone
      #
            jzone  = div(zone - 1, x_zones) + 1
            izone  = mod(zone - 1, x_zones) + 1
      #
            comm_index = 0
            if (zone_proc_id[iz_west[zone]] == iproc)
               comm_index = comm_index + y_size[jzone]
            end
            if (zone_proc_id[iz_east[zone]] == iproc)
               comm_index = comm_index + y_size[jzone]
            end
            if (zone_proc_id[iz_south[zone]] == iproc)
               comm_index = comm_index + x_size[izone]
            end
            if (zone_proc_id[iz_north[zone]] == iproc)
               comm_index = comm_index + x_size[izone]
            end
      #
            return comm_index
end

function map_zones(num_clusters, num_zones, zone_mapping)
   proc_num_zones = Array{Int64}(undef, num_clusters)
   zone_proc_id = zeros(Int64, num_zones)
   proc_zone_id = Array{Array{Int64}}(undef, num_clusters)  # [zeros(Int64, num_zones) for _ = 1:num_clusters]

   for (c, zs) in zone_mapping
       proc_zone_id[c-1] = zs
       proc_num_zones[c-1] = length(zs)
       for z in zs
           zone_proc_id[z] = c-2
       end
   end

   return proc_num_zones, proc_zone_id, zone_proc_id
end

function map_zones(num_clusters, x_zones, y_zones, num_zones, nx, ny, nz, mz_bload, tot_threads, npb_verbose)
      #
      #  Perform zone-process mapping for load balance
      #
      # ... sort the zones in decending order

      proc_num_zones = Array{Int64}(undef, num_clusters)
      zone_proc_id = zeros(Int64, num_zones)
      proc_zone_id = [zeros(Int64, num_zones) for _ = 1:num_clusters]

      z_order = Array{Int64}(undef, max_zones)
      zone_size = Array{FloatType}(undef, max_zones)
      proc_group_flag = Array{Int64}(undef, num_clusters)

      tot_size = 0.0e0
      for iz = 1:num_zones
         zone_size[iz] = 1.0e0*nx[iz]*ny[iz]*nz[iz]
         z_order[iz] = iz
         tot_size = tot_size + zone_size[iz]
      end
      for iz = 1:num_zones-1
         cur_size = zone_size[z_order[iz]]
         mz = iz
         for z2 = iz+1:num_zones
            if cur_size < zone_size[z_order[z2]]
               cur_size = zone_size[z_order[z2]]
               mz = z2
            end
         end
         if mz != iz
            z2 = z_order[iz]
            z_order[iz] = z_order[mz]
            z_order[mz] = z2
         end
      end
#
      if npb_verbose > 1 
         @printf(stdout, "\n Sorted zones:\n  seq. zone    nx    ny    nz    size\n", )
         for iz = 1:num_zones
            z2 = z_order[iz]
            @printf(stdout, "%5i: %5i %5i %5i %5i %9.0F\n", iz, z2, nx[z2], ny[z2], nz[z2], zone_size[z2])
         end

      end

#   10 format(/' Sorted zones:'/  		             '  seq. zone    nx    ny    nz    size')
#   15 format(i5,':',4(1x,i5),1x,f9.0)
#
# ... balance the load among processes
      for ip = 1:num_clusters
         proc_zone_count[ip] = 0
         proc_zone_size[ip] = 0.0e0
      end

      if (abs(mz_bload) > 1) @goto L110 end
#
# ... try a simple block packing scheme
      np1 = np2 = 0
      n1 = 1
      n2 = num_clusters
      work1 = mod(y_zones, num_clusters)
      while (n1 <= n2)
         if n1*n2 == num_clusters
            work2 = mod(x_zones, n1) + mod(y_zones, n2)
            if work2 <= work1
               work1 = work2
               np1 = n1
               np2 = n2
            end
         end
         n1 = n1 + 1
         n2 = div(num_clusters, n1)
      end
      
# ... if can't find a good solution, fall back to next method
      if (work1 != 0) @goto L110 end
#
# ... amount of work for each block
      work1 = div(x_zones, np1)
      work2 = div(y_zones, np2)
#
      for iz = 1:num_zones
         n1 = div(mod(iz-1, x_zones), work1)
         n2 = div(div(iz-1, x_zones), work2)
         ip = n1 + n2*np1 + 1
#
#  ...   assign the zone to the current processor group
         zone_proc_id[iz] = ip - 1
         proc_zone_size[ip] = proc_zone_size[ip] + zone_size[iz]
         proc_zone_count[ip] = proc_zone_count[ip] + 1
      end
#
         @goto L150
#
# ... use a bin-packing scheme to balance the load among processes
      @label L110
      for iz = 1:num_zones
         zone_proc_id[iz] = -1
      end

      iz = 1
      while (iz <= num_zones)
#
#  ...   the current most empty processor
         np = 1
         cur_size = proc_zone_size[1]
         for ip = 2:num_clusters
            if cur_size > proc_zone_size[ip]
               np = ip
               cur_size = proc_zone_size[ip]
            end
         end
         ip = np - 1
#
#  ...   get a zone that has the largest communication index with
#        the current group and does not worsen the computation balance
         mz = z_order[iz]
         if iz < num_zones
            zone_comm = get_comm_index(mz, ip, x_zones, y_zones, zone_proc_id, x_size, y_size)
            for z2 = iz+1:num_zones
               zone = z_order[z2]

               diff_ratio = (zone_size[z_order[iz]] - zone_size[zone]) / zone_size[z_order[iz]]
               if (diff_ratio > 0.05E0) @goto L120 end

               if zone_proc_id[zone] < 0
                  comm_index = get_comm_index(zone, ip, x_zones, y_zones, zone_proc_id, x_size, y_size)
                  if comm_index > zone_comm
                     mz = zone
                     zone_comm = comm_index
                  end
               end
            end
         end

         #
#  ...   assign the zone to the current processor group
         @label L120
         zone_proc_id[mz] = ip
         proc_zone_size[np] = proc_zone_size[np] + zone_size[mz]
         proc_zone_count[np] = proc_zone_count[np] + 1
#
#  ...   skip the previously assigned zones
         while (iz <= num_zones)
            if (zone_proc_id[z_order[iz]] < 0) @goto L130 end
            iz = iz + 1
         end
         @label L130
      end
#
# ... move threads around if needed
      @label L150
      mz = 1
      if (tot_threads <= num_clusters || mz_bload < 1) mz = 0 end
#
      if mz != 0
#
         for ipg = 1:num_clusters
            proc_group_flag[ipg] = 0
         end
#
         ipg = 1
#
# ...    balance load within a processor group
         @label L200
         while (ipg <= num_clusters)
            if (proc_group_flag[ipg] == 0) @goto L210 end
            ipg = ipg + 1
         end
         @label L210
         if (ipg > num_clusters) @goto L300 end
#

         group = proc_group[ipg]
         tot_group_size = 0.0e0
         tot_group_threads = 0
         for ip = ipg:num_clusters
            if proc_group[ip] == group
               proc_group_flag[ip] = 1
               tot_group_size = tot_group_size + proc_zone_size[ip]
               tot_group_threads = tot_group_threads + num_processes[ip]
            end
         end
#
         ave_size = div(tot_group_size, tot_group_threads)
#
#  ...   distribute size evenly among threads
         icur_size = 0
         for ip = 1:num_clusters
            if (proc_group[ip] != group) @goto L220 end
            num_processes[ip] = div(proc_zone_size[ip], ave_size)
            if (num_processes[ip] < 1)
                  num_processes[ip] = 1
            end
            if (max_threads > 0 &&
                  num_processes[ip] > max_threads)
                  num_processes[ip] = max_threads
            end
            icur_size = icur_size + num_processes[ip]
            @label L220
         end
         mz = tot_group_threads - icur_size
#
#  ...   take care of any remainers
         inc = 1
         if (mz < 0) inc = -1 end
         while (mz != 0)
            max_size = 0.0e0
            imx = 0
            for ip = 1:num_clusters
               if (proc_group[ip] != group) @goto L230 end
               if mz > 0
                  cur_size = div(proc_zone_size[ip], num_processes[ip])
                  if cur_size > max_size && (max_threads <= 0||
                        num_processes[ip] < max_threads)
                     max_size = cur_size
                     imx = ip
                  end
               elseif num_processes[ip] > 1
                  cur_size = div(proc_zone_size[ip], num_processes[ip]-1)
                  if max_size == 0 || cur_size < max_size
                     max_size = cur_size
                     imx = ip
                  end
               end
               @label L230
            end
            num_processes[imx] = num_processes[imx] + inc
            mz = mz - inc
         end
#
            @goto L200
      end
#

# ... print the mapping
      @label L300
      if npb_verbose > 0 
         @printf(stdout, "\n Zone-process mapping:\n  proc  nzones  zone_size nthreads size_per_thread\n", )
         for ip = 1:num_clusters
            @printf(stdout, "%6i  %5i  %10.0F  %5i   %10.0F\n", ip-1, proc_zone_count[ip],
                              proc_zone_size[ip], num_processes[ip], div(proc_zone_size[ip],num_processes[ip]))
            for iz = 1:num_zones
               if zone_proc_id[iz] == ip-1
                  @printf(stdout, "   zone  %5i  %9.0F\n", iz, zone_size[iz])
               end
            end
         end
      end

      imx = 1
      max_size = div(proc_zone_size[1], num_processes[1])
      imn = imx
      ave_size = max_size
      for ip = 2:num_clusters
         cur_size = div(proc_zone_size[ip], num_processes[ip])
         if cur_size > max_size
            imx = ip
            max_size = cur_size
         end
         if cur_size < ave_size
            imn = ip
            ave_size = cur_size
         end
      end

      if npb_verbose > 0
         println(stdout, )
         @printf(stdout, " %s: proc = %6i nzones = %5i size = %10.0F nthreads = %5i\n", "Max", imx-1, proc_zone_count[imx],
                     proc_zone_size[imx], num_processes[imx])
         @printf(stdout, " %s: proc = %6i nzones = %5i size = %10.0F nthreads = %5i\n", "Min", imn-1, proc_zone_count[imn],
                     proc_zone_size[imn], num_processes[imn])
      end
      @printf(stdout, "\n Calculated speedup = %9.2F\n\n", tot_size / max_size)

      for clusterid = 0:num_clusters-1
   #
   # ... reorganize list of zones for this process
         zone = 0
         for iz = 1:num_zones
            if zone_proc_id[iz] == clusterid
               zone = zone + 1
               proc_zone_id[clusterid+1][zone] = iz
            end
         end
         proc_num_zones[clusterid+1] = zone
         if zone != proc_zone_count[clusterid+1]
            println(stdout, "Warning: ", clusterid, ": mis-matched zone counts -", zone, proc_zone_count[clusterid+1])
         end
   #
   # ... set number of threads for this process
         group = proc_group[clusterid+1]
         np = 0
         for ip = 1:num_clusters
            if proc_group[ip] == group
               proc_group[ip] = np
               np = np + 1
               num_processes[np] = num_processes[ip]
            end
         end
         ipg = proc_group[clusterid+1]
         if npb_verbose > 1 #&& node == root
            @printf(stdout, " cluster id%6i group%5i group_size%5i group_pid%5i threads%4i\n", clusterid, group, np, ipg, num_processes[ipg+1])
         end 
      end
#
# ... pin-to-node within one process group
#      call smp_pinit_thread(np, ipg, num_processes)
#
      return proc_num_zones, proc_zone_id, zone_proc_id
end
      