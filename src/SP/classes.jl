
@enum CLASS CLASS_UNDEFINED=-1 CLASS_S=0 CLASS_W=1 CLASS_A=2 CLASS_B=3 CLASS_C=4 CLASS_D=5 CLASS_E=6 CLASS_F=7 CLASS_X1=8

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

const sp_class = Dict([ (CLASS_UNDEFINED, Params(100, 0.015e0, 1.0, 4, 4, 24, 24, 6)),
                        (CLASS_S, Params(100, 0.015e0, 1.0, 1, 1, 24, 24, 6 )),
                        (CLASS_W, Params(400, 0.0015e0, 1.0, 4, 4, 64, 64, 8)),
                        (CLASS_A, Params(400, 0.0015e0, 1.0, 4, 4, #=4, 4,=# 128, 128, 16)),
                        (CLASS_B, Params(400, 0.001e0, 1.0, 8, 8, 304, 208, 17)),
                        (CLASS_C, Params(400, 0.00067e0, 1.0, 2, 1, 480, 320, 28)),
                        (CLASS_D, Params(500, 0.0003e0, 1.0, 3, 1#=32, 32=#, 1632, 1216, 34)),
                        (CLASS_E, Params(500, 0.0002e0, 1.0, 2, 2, 4224, 3456, 92)),     # 1.343.029.248
                        (CLASS_X1, Params(500, 0.0002e0, 1.0, 1, 1, 1024, 512, 32 #=64, 64, 4000, 3200, 100=#)),    # 
                        (CLASS_F, Params(500, 0.0001e0, 1.0, 128, 128, 12032, 8960, 250)),
                        ])

export CLASS
