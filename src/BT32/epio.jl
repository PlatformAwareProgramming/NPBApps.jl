#---------------------------------------------------------------------
#---------------------------------------------------------------------

 function setup_btio()

      if node < 10000
          @printf(newfilenm, "%s.%04i\n", filenm, node)
      else
          println(stdout, "error generating file names (> 10000 nodes)")
          exit(1)
      end

      @label L996
      format(a, '.', "i4.4")

      OPEN(unit = 99, file = newfilenm, form = "unformatted",
            status = "unknown")

      for m = 1:5
         xce_sub[m] = 0.0e0
      end

      idump_sub = 0

      return nothing
end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

function output_timestep()

      for cio = 1:ncells
          println(99, ((((u[aio, ix, jio, kio, cio], aio = 1, 5),
                   ix = 0, cell_size[z][1, cio]-1),
                   jio = 0, cell_size[z][2, cio]-1),
                   kio = 0, cell_size[z][3, cio]-1))
      end

      idump_sub = idump_sub + 1
      if rd_interval > 0
         if idump_sub >= rd_interval

            rewind(64)
            acc_sub_norms(idump+1)

            rewind(64)
            idump_sub = 0
         end
      end

      return nothing
end

#---------------------------------------------------------------------
#---------------------------------------------------------------------

function acc_sub_norms(idump_cur)

      ichunk = idump_cur - idump_sub + 1
      for ii = 0:idump_sub-1
        for cio = 1:ncells
          READ(99, ((((u[m, ix, jio, kio, cio], m = 1, 5),
                   ix = 0, cell_size[z][1, cio]-1),
                   jio = 0, cell_size[z][2, cio]-1),
                   kio = 0, cell_size[z][3, cio]-1))
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

      CLOSE(unit = 99)

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

#      character(128) newfilenm
#      integer m

      if (rd_interval > 0) @goto L20 end

      if node < 10000
          @printf(newfilenm, "%s.%04i\n", filenm, node)
      else
          println(stdout, "error generating file names (> 10000 nodes)")
          exit(1)
      end

      @label L996
      format(a, '.', "i4.4")

      OPEN(unit = 99, file = newfilenm,
            form = "unformatted", action = "read")

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
