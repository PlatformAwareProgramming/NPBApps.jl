#---------------------------------------------------------------------
#---------------------------------------------------------------------

function compute_buffer_size(no_nodes, ncells, cell_coord, cell_size, dim)

    if (no_nodes > 1 && ncells > 1)

#---------------------------------------------------------------------
#      compute the actual sizes of the buffers; note that there is 
#      always one cell face that doesn't need buffer space, because it 
#      is at the boundary of the grid
#---------------------------------------------------------------------

        west_size = 0
        east_size = 0

        for c = 1:ncells
           face_size = cell_size[2, c] * cell_size[3, c] * dim * 2
           if (cell_coord[1, c] != 1) 
             west_size = west_size + face_size 
           end
           if (cell_coord[1, c] != ncells) 
             east_size = east_size + face_size
           end
        end

        north_size = 0
        south_size = 0
        for c = 1:ncells
           face_size = cell_size[1, c]*cell_size[3, c] * dim * 2
           if (cell_coord[2, c] != 1) 
             south_size = south_size + face_size 
           end
           if (cell_coord[2, c] != ncells) 
             north_size = north_size + face_size
           end
        end

        top_size = 0
        bottom_size = 0
        for c = 1:ncells
           face_size = cell_size[1, c] * cell_size[2, c] * dim * 2
           if (cell_coord[3, c] != 1) 
             bottom_size = bottom_size + face_size
           end
           if (cell_coord[3, c] != ncells) 
             top_size = top_size + face_size
           end
        end

        start_send_west   = 1
        start_send_east   = start_send_west   + west_size
        start_send_south  = start_send_east   + east_size
        start_send_north  = start_send_south  + south_size
        start_send_bottom = start_send_north  + north_size
        start_send_top    = start_send_bottom + bottom_size
        start_recv_west   = 1
        start_recv_east   = start_recv_west   + west_size
        start_recv_south  = start_recv_east   + east_size
        start_recv_north  = start_recv_south  + south_size
        start_recv_bottom = start_recv_north  + north_size
        start_recv_top    = start_recv_bottom + bottom_size
        ss = SA[start_send_east::Int start_send_west::Int start_send_north::Int start_send_south::Int start_send_top::Int start_send_bottom::Int]
        sr = SA[start_recv_east::Int start_recv_west::Int start_recv_north::Int start_recv_south::Int start_recv_top::Int start_recv_bottom::Int]
        b_size = SA[east_size::Int west_size::Int north_size::Int south_size::Int top_size::Int bottom_size::Int]
    else
        ss = SA[0]
        sr = SA[0]
        b_size = SA[0]
    end

    return ss, sr, b_size
       
end

