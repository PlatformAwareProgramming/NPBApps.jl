#---------------------------------------------------------------------
#---------------------------------------------------------------------



function adi(z, no_nodes,
              ncells::Val{nc},
              slice,
              cell_size,
              cell_start,
              cell_end,
              cell_coord,
              cell_low,
              cell_high,
              u,
              rhs,
              lhs,
              rho_i,
              us,
              cv,
              rhoq,
              rhon,
              rhos,
              vs,
              ws,
              square,
              qs,
              ainv,
              speed,
              forcing,
              dt,
              tx2,
              ty2,
              tz2,
              c1,
              c2,
              c1c2,
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
              dssp,
              con43,
              dx2, dx5, c3c4, c1c5,  c2dtty1, dxmax, dx1, dttx1, dttx2,
              comz5, comz4, comz1, comz6,
              dy3, dy5, dy1, dymax, dtty1, dtty2, 
              dz4, dz5, dz1, dttz2, dttz1, dzmax,
              bt,
              successor,
              predecessor,
              in_buffer,
              out_buffer,
              comm_solve,
              comm_rhs,
              ss,
              sr,
              b_size,
              s,
              requests
       ) where nc

         copy_faces(true, z,
                     no_nodes,
                     ncells,
                     successor,
                     predecessor, 
                     cell_size,
                     cell_start,
                     cell_end,
                     cell_coord,
                     cell_low,
                     cell_high,
                     u,
                     rhs,
                     rho_i,
                     us,
                     vs,
                     ws,
                     square,
                     qs,
                     ainv,
                     speed,
                     forcing,
                     dt,
                     tx2,
                     ty2,
                     tz2,
                     c1,
                     c2,
                     c1c2,
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
                     dssp,
                     con43,
                     in_buffer,
                     out_buffer,
                     requests,
                     timeron,
                     comm_rhs,
                     ss,
                     sr,
                     b_size,
                     )

       txinvr(ncells,
              cell_size,
              cell_start,
              cell_end,
              rho_i,
              us,
              vs,
              ws,
              qs,
              speed,
              ainv,
              rhs,
              c2,
              bt)

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
              rho_i,
              us,
              rhon,
              cv,
              speed,
              dx2, dx5, con43, c3c4, c1c5, c2dttx1, dxmax, dx1, dttx1, dttx2,
              comz5, comz4, comz1, comz6,
              in_buffer,
              out_buffer,
              comm_solve,
              requests,
              s,
              timeron
              )

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
                 rho_i,
                 vs,
                 rhoq,
                 cv,
                 speed,
                 con43, c3c4, c1c5, c2dtty1,
                 dy3, dy5, dy1, dymax, dtty1, dtty2, 
                 comz5, comz4, comz1, comz6,
                 in_buffer,
                 out_buffer,
                 comm_solve,
                 requests,
                  s,
                  timeron
                 )


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
              u,
              us,
              vs,
              qs,
              ainv,
              rho_i,
              ws,
              rhos,
              cv,
              speed,
              dz4, dz5, dz1, con43, c3c4, c1c5, dttz2, dttz1, dzmax,
              comz5, comz4, comz1, comz6, c2dttz1,
              in_buffer,
              out_buffer,
              comm_solve,
              requests,
              s,
              timeron
              ) 

         add(ncells,
             cell_size,
             cell_start,
             cell_end,
             u,
             rhs)

       return nothing
end

