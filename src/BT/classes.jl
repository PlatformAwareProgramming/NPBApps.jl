
@enum CLASS CLASS_UNDEFINED=-1 CLASS_S=0 CLASS_W=1 CLASS_A=2 CLASS_B=3 CLASS_C=4 CLASS_D=5 CLASS_E=6 CLASS_F=7 

struct Params
    problem_size
    niter
    dt
end

const bt_class = Dict([ (CLASS_UNDEFINED, Params(12, 60, 0.010e0)),
                        (CLASS_S, Params(12, 60, 0.010e0)),
                        (CLASS_W, Params(24, 200, 0.0008e0)),
                        (CLASS_A, Params(64, 200, 0.0008e0)),
                        (CLASS_B, Params(102, 200, 0.0003e0)),
                        (CLASS_C, Params(162 ,200, 0.0001e0)),
                        (CLASS_D, Params(408, 250, 0.00002e0)),
                        (CLASS_E, Params(1020, 250, 0.4e-5)),
                        (CLASS_F, Params(2560, 250, 0.6e-6)),
                        ])

export CLASS