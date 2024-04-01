#---------------------------------------------------------------------
#---------------------------------------------------------------------


const c1 = 1.4e0
const c2 = 0.4e0
const c3 = 0.1e0
const c4 = 1.0e0
const c5 = 1.4e0

const bt = sqrt(0.5e0)

const c1c2 = c1 * c2
const c1c5 = c1 * c5
const c3c4 = c3 * c4
const c1345 = c1c5 * c3c4

const conz1 = (1.0e0-c1c5)

const dx1 = 0.75e0
const dx2 = 0.75e0
const dx3 = 0.75e0
const dx4 = 0.75e0
const dx5 = 0.75e0

const dy1 = 0.75e0
const dy2 = 0.75e0
const dy3 = 0.75e0
const dy4 = 0.75e0
const dy5 = 0.75e0

const dz1 = 1.0e0
const dz2 = 1.0e0
const dz3 = 1.0e0
const dz4 = 1.0e0
const dz5 = 1.0e0

const dxmax = max(dx3, dx4)
const dymax = max(dy2, dy4)
const dzmax = max(dz2, dz3)

const dssp = 0.25e0 * max(dx1, max(dy1, dz1) )

const c4dssp = 4.0e0 * dssp
const c5dssp = 5.0e0 * dssp

const c2iv  = 2.5e0
const con43 = 4.0e0/3.0e0
const con16 = 1.0e0/6.0e0

const ce = SA_F64[2.0e0 0.0e0 0.0e0 4.0e0 5.0e0 3.0e0 0.5e0 0.02e0 0.01e0 0.03e0 0.5e0 0.4e0 0.3e0;
                  1.0e0 0.0e0 0.0e0 0.0e0 1.0e0 2.0e0 3.0e0 0.01e0 0.03e0 0.02e0 0.4e0 0.3e0 0.5e0;
                  2.0e0 2.0e0 0.0e0 0.0e0 0.0e0 2.0e0 3.0e0 0.04e0 0.03e0 0.05e0 0.3e0 0.5e0 0.4e0;
                  2.0e0 2.0e0 0.0e0 0.0e0 0.0e0 2.0e0 3.0e0 0.03e0 0.05e0 0.04e0 0.2e0 0.1e0 0.3e0;
                  5.0e0 4.0e0 3.0e0 2.0e0 0.1e0 0.4e0 0.3e0 0.05e0 0.04e0 0.03e0 0.1e0 0.3e0 0.2e0]

function set_constants(dt, grid_points)

      dnxm1 = 1.0e0 / float(grid_points[1]-1)
      dnym1 = 1.0e0 / float(grid_points[2]-1)
      dnzm1 = 1.0e0 / float(grid_points[3]-1)

      tx1 = 1.0e0 / (dnxm1 * dnxm1)
      tx2 = 1.0e0 / (2.0e0 * dnxm1)
      tx3 = 1.0e0 / dnxm1

      ty1 = 1.0e0 / (dnym1 * dnym1)
      ty2 = 1.0e0 / (2.0e0 * dnym1)
      ty3 = 1.0e0 / dnym1

      tz1 = 1.0e0 / (dnzm1 * dnzm1)
      tz2 = 1.0e0 / (2.0e0 * dnzm1)
      tz3 = 1.0e0 / dnzm1


      dttx1 = dt*tx1
      dttx2 = dt*tx2
      dtty1 = dt*ty1
      dtty2 = dt*ty2
      dttz1 = dt*tz1
      dttz2 = dt*tz2

      c2dttx1 = 2.0e0*dttx1
      c2dtty1 = 2.0e0*dtty1
      c2dttz1 = 2.0e0*dttz1

      dtdssp = dt*dssp

      comz1  = dtdssp
      comz4  = 4.0e0*dtdssp
      comz5  = 5.0e0*dtdssp
      comz6  = 6.0e0*dtdssp

      c3c4tx3 = c3c4*tx3
      c3c4ty3 = c3c4*ty3
      c3c4tz3 = c3c4*tz3

      dx1tx1 = dx1*tx1
      dx2tx1 = dx2*tx1
      dx3tx1 = dx3*tx1
      dx4tx1 = dx4*tx1
      dx5tx1 = dx5*tx1

      dy1ty1 = dy1*ty1
      dy2ty1 = dy2*ty1
      dy3ty1 = dy3*ty1
      dy4ty1 = dy4*ty1
      dy5ty1 = dy5*ty1

      dz1tz1 = dz1*tz1
      dz2tz1 = dz2*tz1
      dz3tz1 = dz3*tz1
      dz4tz1 = dz4*tz1
      dz5tz1 = dz5*tz1

 
      xxcon1 = c3c4tx3*con43*tx3
      xxcon2 = c3c4tx3*tx3
      xxcon3 = c3c4tx3*conz1*tx3
      xxcon4 = c3c4tx3*con16*tx3
      xxcon5 = c3c4tx3*c1c5*tx3

      yycon1 = c3c4ty3*con43*ty3
      yycon2 = c3c4ty3*ty3
      yycon3 = c3c4ty3*conz1*ty3
      yycon4 = c3c4ty3*con16*ty3
      yycon5 = c3c4ty3*c1c5*ty3

      zzcon1 = c3c4tz3*con43*tz3
      zzcon2 = c3c4tz3*tz3
      zzcon3 = c3c4tz3*conz1*tz3
      zzcon4 = c3c4tz3*con16*tz3
      zzcon5 = c3c4tz3*c1c5*tz3

      return dnxm1, dnym1, dnzm1 , 
             tx1, tx2, tx3, 
             ty1, ty2, ty3, 
             tz1, tz2, tz3, 
             dttx1 , dttx2, 
             dtty1, dtty2, 
             dttz1, dttz2, 
             c2dttx1, c2dtty1, c2dttz1, dtdssp, 
             comz1, comz4, comz5, comz6, 
             c3c4tx3, c3c4ty3, c3c4tz3, 
             dx1tx1, dx2tx1, dx3tx1, dx4tx1, dx5tx1,
             dy1ty1, dy2ty1, dy3ty1, dy4ty1, dy5ty1, 
             dz1tz1, dz2tz1, dz3tz1, dz4tz1, dz5tz1 , 
             xxcon1, xxcon2, xxcon3, xxcon4, xxcon5,
             yycon1, yycon2, yycon3, yycon4, yycon5,
             zzcon1, zzcon2, zzcon3, zzcon4, zzcon5

end
