#---------------------------------------------------------------------
#---------------------------------------------------------------------

const send_id = Ref{MPI.Request}(MPI.REQUEST_NULL)
const recv_id = Ref{MPI.Request}(MPI.REQUEST_NULL)

function adi(iz, ss, 
            sr, 
            b_size,
            MAX_CELL_DIM,
            IMAX,
            JMAX,
            KMAX,
            cell_coord,
            cell_size,
            cell_start,
            cell_end,
            slice,
            forcing,           
            u,
            rhs,
            lhsc,
            backsub_info,
            in_buffer,
            out_buffer,
            fjac,
            njac,
            lhsa,
            lhsb,
            us,
            vs,
            ws,
            qs,
            rho_i,
            square,
            dt,
            timeron,
            ncells,
            tx1,
            tx2,
            ty1,
            ty2,
            tz1,
            tz2,
            dx1tx1,
            dx2tx1,
            dx3tx1,
            dx4tx1,
            dx5tx1,
            dy1ty1,
            dy2ty1,
            dy3ty1,
            dy4ty1,
            dy5ty1,
            dz1tz1,
            dz2tz1,
            dz3tz1,
            dz4tz1,
            dz5tz1,
            xxcon2,
            xxcon3,
            xxcon4,
            xxcon5,
            yycon2,
            yycon3,
            yycon4,
            yycon5,
            zzcon2,
            zzcon3,
            zzcon4,
            zzcon5,
            no_nodes, 
            comm_solve,
            comm_rhs,
            predecessor,
            successor,
            utmp,
            requests,
            )

      copy_faces( ss, 
                  sr, 
                  b_size,
                  cell_coord,
                  cell_size,
                  cell_start,
                  cell_end,
                  forcing,        
                  u,
                  rhs,
                  in_buffer,
                  out_buffer,
                  us,
                  vs,
                  ws,
                  qs,
                  rho_i,
                  square,
                  timeron,
                  dt,
                  ncells,
                  tx2,
                  ty2,
                  tz2,
                  dx1tx1,
                  dx2tx1,
                  dx3tx1,
                  dx4tx1,
                  dx5tx1,
                  dy1ty1,
                  dy2ty1,
                  dy3ty1,
                  dy4ty1,
                  dy5ty1,
                  dz1tz1,
                  dz2tz1,
                  dz3tz1,
                  dz4tz1,
                  dz5tz1,
                  xxcon2,
                  xxcon3,
                  xxcon4,
                  xxcon5,
                  yycon2,
                  yycon3,
                  yycon4,
                  yycon5,
                  zzcon2,
                  zzcon3,
                  zzcon4,
                  zzcon5,
                  no_nodes, 
                  comm_rhs,
                  predecessor,
                  successor,
                  requests,
            )


        if iz == 1 
            write_u(iz)
         end

      x_solve(
            MAX_CELL_DIM,
            JMAX,
            KMAX,
            cell_coord,
            cell_size,
            cell_start,
            cell_end,
            slice,
            u,
            rhs,
            lhsc,
            backsub_info,
            in_buffer,
            out_buffer,
            fjac,
            njac,
            lhsa,
            lhsb,
            qs,
            rho_i,
            dt,
            timeron,
            ncells,
            tx1,
            tx2,
            comm_solve,
            predecessor,
            successor
        ) #maxdepth=3 modules=[BT]
   
        y_solve(
            MAX_CELL_DIM,
            IMAX,
            JMAX,
            KMAX,
            cell_coord,
            cell_size,
            cell_start,
            cell_end,
            slice,
            u,
            rhs,
            lhsc,
            backsub_info,
            in_buffer,
            out_buffer,
            fjac,
            njac,
            lhsa,
            lhsb,
            qs,
            dt,
            timeron,
            ncells,
            ty1,
            ty2,
            comm_solve,     
            predecessor,
            successor,
            utmp
      )

      z_solve(
            MAX_CELL_DIM,
            IMAX,
            JMAX,
            cell_coord,
            cell_size,
            cell_start,
            cell_end,
            slice,     
            u,
            rhs,
            lhsc,
            backsub_info,
            in_buffer,
            out_buffer,
            fjac,
            njac,
            lhsa,
            lhsb,
            qs,
            dt,
            timeron,
            ncells,
            tz1,
            tz2,
            comm_solve,
            predecessor,
            successor,
            utmp
     )

      add(cell_size,
            cell_start,
            cell_end,         
            u,
            rhs,
            ncells,   
            )

      return nothing
end

