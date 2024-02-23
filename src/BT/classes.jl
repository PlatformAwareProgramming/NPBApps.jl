@enum CLASS CLASS_UNDEFINED=-1 CLASS_S=0 CLASS_W=1 CLASS_A=2 CLASS_B=3 CLASS_C=4 CLASS_D=5 CLASS_E=6 CLASS_F=7 

struct Params
    niter
    dt
    ratio
    x_zones
    y_zones
    gx_size
    gy_size
    gz_size
end

const sp_class = Dict([ (CLASS_UNDEFINED, Params(100, 0.010e0, 3.0, 2, 2, 24, 24, 6)),
                        (CLASS_S, Params(60, 0.010e0, 3.0, 2, 2, 24, 24, 6 )),
                        (CLASS_W, Params(200, 0.0008e0, 4.5, 4, 4, 64, 64, 8)),
                        (CLASS_A, Params(200, 0.0008e0, 4.5, 4, 4, 128, 128, 16)),
                        (CLASS_B, Params(200, 0.0003e0, 4.5, 8, 8, 304, 208, 17)),
                        (CLASS_C, Params(200, 0.0001e0, 4.5, 16, 16, 480, 320, 28)),
                        (CLASS_D, Params(250, 0.00002e0, 4.5, 32, 32, 1632, 1216, 34)),
                        (CLASS_E, Params(250, 0.000004e0, 4.5, 64, 64, 4224, 3456, 92)),
                        (CLASS_F, Params(250, 0.000001e0, 4.5, 128, 128, 12032, 8960, 250)),
                        ])

export CLASS