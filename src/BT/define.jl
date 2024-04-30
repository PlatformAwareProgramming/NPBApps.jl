#---------------------------------------------------------------------
#---------------------------------------------------------------------

function compute_buffer_size_initial(proc_num_zones)

  global west_size   = Array{Int}(undef, proc_num_zones)
  global east_size   = Array{Int}(undef, proc_num_zones)
  global north_size  = Array{Int}(undef, proc_num_zones)
  global south_size  = Array{Int}(undef, proc_num_zones)
  global top_size    = Array{Int}(undef, proc_num_zones)
  global bottom_size = Array{Int}(undef, proc_num_zones)
  
  global start_send_west   = Array{Int}(undef, proc_num_zones)
  global start_send_east   = Array{Int}(undef, proc_num_zones)
  global start_send_south  = Array{Int}(undef, proc_num_zones)
  global start_send_north  = Array{Int}(undef, proc_num_zones)
  global start_send_bottom = Array{Int}(undef, proc_num_zones)
  global start_send_top    = Array{Int}(undef, proc_num_zones)
  global start_recv_west   = Array{Int}(undef, proc_num_zones)
  global start_recv_east   = Array{Int}(undef, proc_num_zones)
  global start_recv_south  = Array{Int}(undef, proc_num_zones)
  global start_recv_north  = Array{Int}(undef, proc_num_zones)
  global start_recv_bottom = Array{Int}(undef, proc_num_zones)
  global start_recv_top    = Array{Int}(undef, proc_num_zones)

end

function compute_buffer_size(z, dim)

#---------------------------------------------------------------------
#      compute the actual sizes of the buffers; note that there is 
#      always one cell face that doesn't need buffer space, because it 
#      is at the boundary of the grid
#---------------------------------------------------------------------

   west_size[z] = 0
   east_size[z] = 0

   for c = 1:ncells
      face_size = cell_size[z][2, c] * cell_size[z][3, c] * dim * 2
      west_size[z] = west_size[z] + face_size 
      east_size[z] = east_size[z] + face_size
   end

   north_size[z] = 0
   south_size[z] = 0
   for c = 1:ncells
      face_size = cell_size[z][1, c]*cell_size[z][3, c] * dim * 2
      south_size[z] = south_size[z] + face_size 
      north_size[z] = north_size[z] + face_size
   end

   top_size[z] = 0
   bottom_size[z] = 0
   for c = 1:ncells
      face_size = cell_size[z][1, c] * cell_size[z][2, c] * dim * 2
      bottom_size[z] = bottom_size[z] + face_size
      top_size[z] = top_size[z] + face_size
   end

   start_send_west[z]   = 1
   start_send_east[z]   = start_send_west[z]   + west_size[z]
   start_send_south[z]  = start_send_east[z]   + east_size[z]
   start_send_north[z]  = start_send_south[z]  + south_size[z]
   start_send_bottom[z] = start_send_north[z]  + north_size[z]
   start_send_top[z]    = start_send_bottom[z] + bottom_size[z]
   start_recv_west[z]   = 1
   start_recv_east[z]   = start_recv_west[z]   + west_size[z]
   start_recv_south[z]  = start_recv_east[z]   + east_size[z]
   start_recv_north[z]  = start_recv_south[z]  + south_size[z]
   start_recv_bottom[z] = start_recv_north[z]  + north_size[z]
   start_recv_top[z]    = start_recv_bottom[z] + bottom_size[z]

return nothing
   
end

