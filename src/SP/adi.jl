#---------------------------------------------------------------------
#---------------------------------------------------------------------

function adi()

       copy_faces()
         
       txinvr()

       x_solve(ncells,
       successor,
       predecessor,
       slice,
       cell_size,
       cell_start,
       cell_end,
       cell_coord,
       lhs,
       rhs,
       in_buffer,
       out_buffer,
       comm_solve)

       y_solve(ncells,
       successor,
       predecessor,
       slice,
       cell_size,
       cell_start,
       cell_end,
       cell_coord,
       lhs,
       rhs,
       in_buffer,
       out_buffer,
       comm_solve)

       z_solve(ncells,
       successor,
       predecessor,
       slice,
       cell_size,
       cell_start,
       cell_end,
       cell_coord,
       lhs,
       rhs,
       in_buffer,
       out_buffer,
       comm_solve)

       add()

       return nothing
end

