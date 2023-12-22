
#---------------------------------------------------------------------
#---------------------------------------------------------------------

function setup_btio()

      iseek = 0

      if node == root
          MPI_File_delete(filenm, MPI_INFO_NULL, ierr)
      end

      MPI_Barrier(comm_solve, ierr)

      MPI_File_open(comm_solve,
                filenm,
                MPI_MODE_RDWR + MPI_MODE_CREATE,
                MPI_INFO_NULL,
                fp,
                ierr)

      MPI_File_set_view(fp,
                iseek, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION,
                "native", MPI_INFO_NULL, ierr)

      if ierr != MPI_SUCCESS
          println(stdout, "Error opening file")
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

#      use bt_data
#      use mpinpb

#      implicit none

#      integer COUNT, jio, kio, cio, aio
#      integer ierr
#      integer mstatus[MPI_STATUS_SIZE]

      for cio = 1:ncells
          for kio = 0:cell_size[3, cio]-1
              for jio = 0:cell_size[2, cio]-1
                  iseek = (cell_low[3, cio]+kio) +
                         PROBLEM_SIZE*idump_sub
                  iseek = (cell_low[2, cio]+jio) +
                         PROBLEM_SIZE*iseek
                  iseek = 5*(cell_low[1, cio] +
                         PROBLEM_SIZE*iseek)

                  COUNT = 5*cell_size[1, cio]

                  MPI_File_write_at(fp, iseek,
                        u[1, 0, jio, kio, cio],
                        COUNT, MPI_DOUBLE_PRECISION,
                        mstatus, ierr)

                  if ierr != MPI_SUCCESS
                      println(stdout, "Error writing to file")
                      exit(1)
                  end
              end
          end
      end

      idump_sub = idump_sub + 1
      if rd_interval > 0
         if idump_sub >= rd_interval

            acc_sub_norms(idump+1)

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

#      integer COUNT, jio, kio, cio, ii, m, ichunk
#      integer ierr
#      integer mstatus[MPI_STATUS_SIZE]
#      DOUBLEPRECISION xce_single[5]

      ichunk = idump_cur - idump_sub + 1
      for ii = 0:idump_sub-1
        for cio = 1:ncells
          for kio = 0:cell_size[3, cio]-1
              for jio = 0:cell_size[2, cio]-1
                  iseek = (cell_low[3, cio]+kio) +
                         PROBLEM_SIZE*ii
                  iseek = (cell_low[2, cio]+jio) +
                         PROBLEM_SIZE*iseek
                  iseek = 5*(cell_low[1, cio] +
                         PROBLEM_SIZE*iseek)

                  COUNT = 5*cell_size[1, cio]

                  MPI_File_read_at(fp, iseek,
                        u[1, 0, jio, kio, cio],
                        COUNT, MPI_DOUBLE_PRECISION,
                        mstatus, ierr)

                  if ierr != MPI_SUCCESS
                      println(stdout, "Error reading back file")
                      MPI_File_close(fp, ierr)
                      exit(1)
                  end
              end
          end
        end

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
      MPI_File_set_view(fp,
                iseek, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION,
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

