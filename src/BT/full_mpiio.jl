


#---------------------------------------------------------------------
#---------------------------------------------------------------------

 function setup_btio()

       mpi_bcast(collbuf_nodes, 1, MPI_INTEGER,
                      root, comm_setup, ierr)

       mpi_bcast(collbuf_size, 1, MPI_INTEGER,
                      root, comm_setup, ierr)

       if collbuf_nodes == 0
          info = MPI_INFO_NULL
       else
          println(cb_nodes, collbuf_nodes)
          println(cb_size, collbuf_size)
          MPI_Info_create(info, ierr)
          MPI_Info_set(info, "cb_nodes", cb_nodes, ierr)
          MPI_Info_set(info, "cb_buffer_size", cb_size, ierr)
          MPI_Info_set(info, "collective_buffering", "true", ierr)
       end

       MPI_Type_contiguous(5, MPI_DOUBLE_PRECISION,
                                element, ierr)
       MPI_Type_commit(element, ierr)
       MPI_Type_extent(element, eltext, ierr)

       for c = 1:ncells
#
# Outer array dimensions ar same for every cell
#
           sizes[1] = IMAX+4
           sizes[2] = JMAX+4
           sizes[3] = KMAX+4
#
# 4th dimension is cell number, total of maxcells cells
#
           sizes[4] = maxcells
#
# Internal dimensions of cells can differ slightly between cells
#
           subsizes[1] = cell_size[1, c]
           subsizes[2] = cell_size[2, c]
           subsizes[3] = cell_size[3, c]
#
# Cell is 4th dimension, 1 cell per cell type to handle varying 
# cell sub-array sizes
#
           subsizes[4] = 1

#
# type constructors use 0-based start addresses
#
           starts[1] = 2
           starts[2] = 2
           starts[3] = 2
           starts[4] = c-1

# 
# Create buftype for a cell
#
           MPI_Type_create_subarray(4, sizes, subsizes,
                starts, MPI_ORDER_FORTRAN, element,
                cell_btype[c], ierr)
#
# block length and displacement for joining cells - 
# 1 cell buftype per block, cell buftypes have own displacment
# generated from cell number (4th array dimension)
#
           cell_blength[c] = 1
           cell_disp[c] = 0

       end
#
# Create combined buftype for all cells
#
       MPI_Type_struct(ncells, cell_blength, cell_disp,
                  cell_btype, combined_btype, ierr)
       MPI_Type_commit(combined_btype, ierr)

       for c = 1:ncells
#
# Entire array size
#
           sizes[1] = PROBLEM_SIZE
           sizes[2] = PROBLEM_SIZE
           sizes[3] = PROBLEM_SIZE

#
# Size of c'th cell
#
           subsizes[1] = cell_size[1, c]
           subsizes[2] = cell_size[2, c]
           subsizes[3] = cell_size[3, c]

#
# Starting point in full array of c'th cell
#
           starts[1] = cell_low[1, c]
           starts[2] = cell_low[2, c]
           starts[3] = cell_low[3, c]

           MPI_Type_create_subarray(3, sizes, subsizes,
                starts, MPI_ORDER_FORTRAN,
                element, cell_ftype[c], ierr)
           cell_blength[c] = 1
           cell_disp[c] = 0
       end

       MPI_Type_struct(ncells, cell_blength, cell_disp,
                  cell_ftype, combined_ftype, ierr)
       MPI_Type_commit(combined_ftype, ierr)

       iseek = 0
       if node == root
          MPI_File_delete(filenm, MPI_INFO_NULL, ierr)
       end


      MPI_Barrier(comm_solve, ierr)

       MPI_File_open(comm_solve,
                filenm,
                MPI_MODE_RDWR+MPI_MODE_CREATE,
                MPI_INFO_NULL, fp, ierr)

       if ierr != MPI_SUCCESS
                println(stdout, "Error opening file")
                exit(1)
       end

        MPI_File_set_view(fp, iseek, element,
                combined_ftype, "native", info, ierr)

       if ierr != MPI_SUCCESS
                println(stdout, "Error setting file view")
                exit(1)
       end

      for m = 1:5
         xce_sub[m] = 0.0e0
      end

      idump_sub = 0


      return nothing
end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

 function output_timestep()

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#      use bt_data
#      use mpinpb

#      implicit none

#      integer mstatus[MPI_STATUS_SIZE]
#      integer ierr

      MPI_File_write_at_all(fp, iseek, u,
                                 1, combined_btype, mstatus, ierr)
      if ierr != MPI_SUCCESS
          println(stdout, "Error writing to file")
          exit(1)
      end

      MPI_Type_size(combined_btype, iosize, ierr)
      iseek = iseek + iosize/eltext

      idump_sub = idump_sub + 1
      if rd_interval > 0
         if idump_sub >= rd_interval

            iseek = 0
            acc_sub_norms(idump+1)

            iseek = 0
            idump_sub = 0
         end
      end

      return nothing
end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

 function acc_sub_norms(idump_cur)

#      use bt_data
#      use mpinpb

#      implicit none

#      integer idump_cur

#      integer ii, m, ichunk
#      integer ierr
#      integer mstatus[MPI_STATUS_SIZE]
#      DOUBLEPRECISION xce_single[5]

      ichunk = idump_cur - idump_sub + 1
      for ii = 0:idump_sub-1

        MPI_File_read_at_all(fp, iseek, u,
                                 1, combined_btype, mstatus, ierr)
        if ierr != MPI_SUCCESS
           println(stdout, "Error reading back file")
           MPI_File_close(fp, ierr)
           exit(1)
        end

        MPI_Type_size(combined_btype, iosize, ierr)
        iseek = iseek + iosize/eltext

        if (node == root) println(stdout, "Reading data set ", ii+ichunk) end

        error_norm(xce_single)
        for m = 1:5
           xce_sub[m] = xce_sub[m] + xce_single[m]
        end
      end

      return nothing
end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

 function btio_cleanup()

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#      use bt_data
#      use mpinpb

#      implicit none

#      integer ierr

      MPI_File_close(fp, ierr)

      return nothing
end

#---------------------------------------------------------------------
#---------------------------------------------------------------------


 function accumulate_norms(xce_acc)

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#      use bt_data
#      use mpinpb

#      implicit none

#      DOUBLEPRECISION xce_acc[5]
#      integer m, ierr

      if (rd_interval > 0) @goto L20 end

      MPI_File_open(comm_solve,
                filenm,
                MPI_MODE_RDONLY,
                MPI_INFO_NULL,
                fp,
                ierr)

      iseek = 0
      MPI_File_set_view(fp, iseek, element, combined_ftype,
                "native", MPI_INFO_NULL, ierr)

#     clear the last time step

      clear_timestep

#     read back the time steps and accumulate norms

      acc_sub_norms(idump)

      MPI_File_close(fp, ierr)

      @label L20
      for m = 1:5
         xce_acc[m] = xce_sub[m] / float(idump)
      end

      return nothing
end

