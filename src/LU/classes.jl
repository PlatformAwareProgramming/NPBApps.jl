@enum CLASS CLASS_UNDEFINED=-1 CLASS_S=0 CLASS_W=1 CLASS_A=2 CLASS_B=3 CLASS_C=4 CLASS_D=5 CLASS_E=6 CLASS_F=7 

struct Params
    isiz01
    isiz02
    isiz03
    itmax
    inorm
    dt
end

const bt_class = Dict([ (CLASS_UNDEFINED, Params(-1, -1, -1, -1, -1, 0.0e0)),
                        (CLASS_S, Params(12, 12, 12, 50, 50, 0.5e0)),
                        (CLASS_W, Params(33, 33, 33, 300, 300, 1.5e-3)),
                        (CLASS_A, Params(64, 64, 64, 250, 250, 2.0e0)),
                        (CLASS_B, Params(102, 102, 102, 250, 250, 2.0e0)),
                        (CLASS_C, Params(162, 162, 162, 250, 250, 2.0e0)),
                        (CLASS_D, Params(408, 408, 408, 300, 300, 1.0e0)),
                        (CLASS_E, Params(1020, 1020, 1020, 300, 300, 0.5e0)),
                        (CLASS_F, Params(2560, 2560, 2560, 300, 300, 0.2e0)),
                        ])

export CLASS