@enum CLASS CLASS_UNDEFINED=-1 CLASS_S=0 CLASS_W=1 CLASS_A=2 CLASS_B=3 CLASS_C=4 CLASS_D=5 CLASS_E=6 CLASS_F=7 

struct Params
    x_zones
    y_zones
    gx_size
    gy_size
    gz_size
    itmax
    inorm
    dt
    ratio
end

const lu_class = Dict([ (CLASS_UNDEFINED, Params(4, 4, 24, 24, 6, 50, 50, 0.5e0, 1.0e0)),
                        (CLASS_S, Params(4, 4, 24, 24, 6, 50, 50, 0.5e0, 1.0e0)),
                        (CLASS_W, Params(4, 4, 64, 64, 8, 300, 300, 1.5e-3, 1.0e0)),
                        (CLASS_A, Params(4, 4, 128, 128, 16, 250, 250, 2.0e0, 1.0e0)),
                        (CLASS_B, Params(4, 4, 304, 208, 17, 250, 250, 2.0e0, 1.0e0)),
                        (CLASS_C, Params(4, 4, 480, 320, 28, 250, 250, 2.0e0, 1.0e0)),
                        (CLASS_D, Params(4, 4, 1632, 1216, 34, 300, 300, 1.0e0, 1.0e0)),
                        (CLASS_E, Params(4, 4, 4224, 3456, 92, 300, 300, 0.5e0, 1.0e0)),
                        (CLASS_F, Params(4, 4, 12032, 8960, 250, 300, 300, 0.2e0, 1.0e0)),
                        ])

export CLASS