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

const  c2iv  = 2.5e0
const  con43 = 4.0e0/3.0e0
const  con16 = 1.0e0/6.0e0

const ce = SA_F64[2.0e0 0.0e0 0.0e0 4.0e0 5.0e0 3.0e0 0.5e0 0.02e0 0.01e0 0.03e0 0.5e0 0.4e0 0.3e0;
                  1.0e0 0.0e0 0.0e0 0.0e0 1.0e0 2.0e0 3.0e0 0.01e0 0.03e0 0.02e0 0.4e0 0.3e0 0.5e0;
                  2.0e0 2.0e0 0.0e0 0.0e0 0.0e0 2.0e0 3.0e0 0.04e0 0.03e0 0.05e0 0.3e0 0.5e0 0.4e0;
                  2.0e0 2.0e0 0.0e0 0.0e0 0.0e0 2.0e0 3.0e0 0.03e0 0.05e0 0.04e0 0.2e0 0.1e0 0.3e0;
                  5.0e0 4.0e0 3.0e0 2.0e0 0.1e0 0.4e0 0.3e0 0.05e0 0.04e0 0.03e0 0.1e0 0.3e0 0.2e0]

#---------------------------------------------------------------------
#---------------------------------------------------------------------

 function set_constants()

       global dnxm1 = 1.0e0 / float(grid_points[1]-1)
       global dnym1 = 1.0e0 / float(grid_points[2]-1)
       global dnzm1 = 1.0e0 / float(grid_points[3]-1)
       
       global tx1 = 1.0e0 / (dnxm1 * dnxm1)
       global tx2 = 1.0e0 / (2.0e0 * dnxm1)
       global tx3 = 1.0e0 / dnxm1
       
       global ty1 = 1.0e0 / (dnym1 * dnym1)
       global ty2 = 1.0e0 / (2.0e0 * dnym1)
       global ty3 = 1.0e0 / dnym1
       
       global tz1 = 1.0e0 / (dnzm1 * dnzm1)
       global tz2 = 1.0e0 / (2.0e0 * dnzm1)
       global tz3 = 1.0e0 / dnzm1

       global dttx1 = dt*tx1
       global dttx2 = dt*tx2
       global dtty1 = dt*ty1
       global dtty2 = dt*ty2
       global dttz1 = dt*tz1
       global dttz2 = dt*tz2
       
       global c2dttx1 = 2.0e0*dttx1
       global c2dtty1 = 2.0e0*dtty1
       global c2dttz1 = 2.0e0*dttz1
       
       global dtdssp = dt*dssp
       
       global comz1  = dtdssp
       global comz4  = 4.0e0*dtdssp
       global comz5  = 5.0e0*dtdssp
       global comz6  = 6.0e0*dtdssp
                            
       global  c3c4tx3 = c3c4*tx3
       global c3c4ty3 = c3c4*ty3
       global  c3c4tz3 = c3c4*tz3
       
       global dx1tx1 = dx1*tx1
       global  dx2tx1 = dx2*tx1
       global dx3tx1 = dx3*tx1
       global  dx4tx1 = dx4*tx1
       global  dx5tx1 = dx5*tx1
       
       global dy1ty1 = dy1*ty1
       global dy2ty1 = dy2*ty1
       global dy3ty1 = dy3*ty1
       global  dy4ty1 = dy4*ty1
       global dy5ty1 = dy5*ty1
       
       global  dz1tz1 = dz1*tz1
       global  dz2tz1 = dz2*tz1
       global  dz3tz1 = dz3*tz1
       global  dz4tz1 = dz4*tz1
       global  dz5tz1 = dz5*tz1
       
       global  xxcon1 = c3c4tx3*con43*tx3
       global  xxcon2 = c3c4tx3*tx3
       global  xxcon3 = c3c4tx3*conz1*tx3
       global  xxcon4 = c3c4tx3*con16*tx3
       global  xxcon5 = c3c4tx3*c1c5*tx3
       
       global  yycon1 = c3c4ty3*con43*ty3
       global  yycon2 = c3c4ty3*ty3
       global  yycon3 = c3c4ty3*conz1*ty3
       global  yycon4 = c3c4ty3*con16*ty3
       global  yycon5 = c3c4ty3*c1c5*ty3
       
       global  zzcon1 = c3c4tz3*con43*tz3
       global  zzcon2 = c3c4tz3*tz3
       global  zzcon3 = c3c4tz3*conz1*tz3
       global  zzcon4 = c3c4tz3*con16*tz3
       global  zzcon5 = c3c4tz3*c1c5*tz3  

#=       global ce = zeros(Float64, 5, 13)

       ce[1, 1]  = 2.0e0
       ce[1, 2]  = 0.0e0
       ce[1, 3]  = 0.0e0
       ce[1, 4]  = 4.0e0
       ce[1, 5]  = 5.0e0
       ce[1, 6]  = 3.0e0
       ce[1, 7]  = 0.5e0
       ce[1, 8]  = 0.02e0
       ce[1, 9]  = 0.01e0
       ce[1, 10] = 0.03e0
       ce[1, 11] = 0.5e0
       ce[1, 12] = 0.4e0
       ce[1, 13] = 0.3e0

       ce[2, 1]  = 1.0e0
       ce[2, 2]  = 0.0e0
       ce[2, 3]  = 0.0e0
       ce[2, 4]  = 0.0e0
       ce[2, 5]  = 1.0e0
       ce[2, 6]  = 2.0e0
       ce[2, 7]  = 3.0e0
       ce[2, 8]  = 0.01e0
       ce[2, 9]  = 0.03e0
       ce[2, 10] = 0.02e0
       ce[2, 11] = 0.4e0
       ce[2, 12] = 0.3e0
       ce[2, 13] = 0.5e0

       ce[3, 1]  = 2.0e0
       ce[3, 2]  = 2.0e0
       ce[3, 3]  = 0.0e0
       ce[3, 4]  = 0.0e0
       ce[3, 5]  = 0.0e0
       ce[3, 6]  = 2.0e0
       ce[3, 7]  = 3.0e0
       ce[3, 8]  = 0.04e0
       ce[3, 9]  = 0.03e0
       ce[3, 10] = 0.05e0
       ce[3, 11] = 0.3e0
       ce[3, 12] = 0.5e0
       ce[3, 13] = 0.4e0

       ce[4, 1]  = 2.0e0
       ce[4, 2]  = 2.0e0
       ce[4, 3]  = 0.0e0
       ce[4, 4]  = 0.0e0
       ce[4, 5]  = 0.0e0
       ce[4, 6]  = 2.0e0
       ce[4, 7]  = 3.0e0
       ce[4, 8]  = 0.03e0
       ce[4, 9]  = 0.05e0
       ce[4, 10] = 0.04e0
       ce[4, 11] = 0.2e0
       ce[4, 12] = 0.1e0
       ce[4, 13] = 0.3e0

       ce[5, 1]  = 5.0e0
       ce[5, 2]  = 4.0e0
       ce[5, 3]  = 3.0e0
       ce[5, 4]  = 2.0e0
       ce[5, 5]  = 0.1e0
       ce[5, 6]  = 0.4e0
       ce[5, 7]  = 0.3e0
       ce[5, 8]  = 0.05e0
       ce[5, 9]  = 0.04e0
       ce[5, 10] = 0.03e0
       ce[5, 11] = 0.1e0
       ce[5, 12] = 0.3e0
       ce[5, 13] = 0.2e0
=#

       return nothing
end
