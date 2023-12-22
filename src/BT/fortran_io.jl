#---------------------------------------------------------------------
#---------------------------------------------------------------------

function setup_btio()

#     determine a proper record_length to use
      if node == root
         frec_sz = fortran_rec_sz
         if frec_sz > 0
            # use the compiled value
            record_length = 40/frec_sz
         else
            # query directly
            inquire(iolength = record_length) d5
         end
         if (record_length < 1) record_length = 40 end
      end

      mpi_bcast(record_length, 1, MPI_INTEGER,
                      root, comm_setup, ierr)

      OPEN(unit = 99, file = filenm,
            form = "unformatted", access = "direct",
            recl = record_length)

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
#      implicit none

#      integer ix, jio, kio, cio

      for cio = 1:ncells
          for kio = 0:cell_size[3, cio]-1
              for jio = 0:cell_size[2, cio]-1
                  iseek = (cell_low[3, cio]+kio) +
                         PROBLEM_SIZE*idump_sub
                  iseek = (cell_low[2, cio]+jio) +
                         PROBLEM_SIZE*iseek
                  iseek = (cell_low[1, cio]) +
                         PROBLEM_SIZE*iseek

                  for ix = 0:cell_size[1, cio]-1
                      println(99, u[1, ix, jio, kio, cio],
                            u[2, ix, jio, kio, cio],
                            u[3, ix, jio, kio, cio],
                            u[4, ix, jio, kio, cio],
                            u[5, ix, jio, kio, cio])
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

#      integer ix, jio, kio, cio, ii, m, ichunk
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
                  iseek = (cell_low[1, cio]) +
                         PROBLEM_SIZE*iseek


                  for ix = 0:cell_size[1, cio]-1
                      READ(99, u[1, ix, jio, kio, cio],
                            u[2, ix, jio, kio, cio],
                            u[3, ix, jio, kio, cio],
                            u[4, ix, jio, kio, cio],
                            u[5, ix, jio, kio, cio])
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

#      implicit none

      CLOSE(unit = 99)

      return nothing
end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

 function accumulate_norms(xce_acc)

#      use bt_data
#      implicit none

#      DOUBLEPRECISION xce_acc[5]
#      integer m

      if (rd_interval > 0) @goto L20 end

      OPEN(unit = 99, file = filenm,
            form = "unformatted", access = "direct",
            recl = record_length, action = "read")

#     clear the last time step

      clear_timestep

#     read back the time steps and accumulate norms

      acc_sub_norms(idump)

      CLOSE(unit = 99)

      @label L20
      for m = 1:5
         xce_acc[m] = xce_sub[m] / float(idump)
      end

      return nothing
end
