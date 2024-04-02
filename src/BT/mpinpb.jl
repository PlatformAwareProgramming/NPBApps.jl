

#      use sp_data, ONLY : max_zones, kind2
#
#
#     zone_proc_id(MZ)     - process id each zone assigned to
#     proc_zone_id(MZ)     - list of zones assigned to this process
#     proc_num_zones       - number of zones assigned to this process
#     proc_zone_count(NP)  - number of zones assigned to each process
#     proc_num_threads(NP) - number of threads assigned to each process
#     proc_group(NP)       - group id each process assigned to
#
#      integer   zone_proc_id[max_zones], proc_zone_id[max_zones],  
#                proc_num_zones
#      integer, allocatable :: proc_zone_count[:],  
#                proc_num_threads[:], proc_group[:]
#      DOUBLEPRECISION, allocatable :: proc_zone_size[:]
#
#      integer   clusterid, root, comm_setup, ierror, dp_type
#      integer   num_threads, max_threads, mz_bload, mz_bload_erank
#      integer   num_clusters, num_clusters2
#      logical   active
#
# ... Two adjustable parameters for MPI communication
#     max_reqs  -- max. number of async message requests
#     MSG_SIZE  -- optimal message size (in words) for communication
#      integer   max_reqs, MSG_SIZE
#                 max_reqs = 32; MSG_SIZE = 400000
#
#      integer   requests[max_reqs], statuses[MPI_STATUS_SIZE,max_reqs]
#
#      integer, allocatable :: pcomm_group[:]
#      integer(kind2), allocatable :: qcomm_size[:]



#---------------------------------------------------------------------
#---------------------------------------------------------------------
#
#  Allocate process-based working arrays
#
#---------------------------------------------------------------------
#---------------------------------------------------------------------

function alloc_proc_space(num_clusters, max_zones)

      num_clusters2 = 1
      # ... proc size that is a power of two and no less than num_clusters
      while (num_clusters2 < num_clusters)
         num_clusters2 = num_clusters2 * 2
      end

      global proc_zone_count = zeros(Int64, num_clusters)
      global num_processes = zeros(Int64, num_clusters)
      global proc_group = zeros(Int64, num_clusters)
      global proc_zone_size = zeros(Int64, num_clusters)
      global pcomm_group = zeros(Int64, num_clusters2)

      return nothing
end
