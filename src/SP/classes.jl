
@enum CLASS CLASS_UNDEFINED=-1 CLASS_S=0 CLASS_W=1 CLASS_A=2 CLASS_B=3 CLASS_C=4 CLASS_D=5 CLASS_E=6 CLASS_F=7 

struct Params
    problem_size
    niter
    dt
end

const sp_class = Dict([ (CLASS_UNDEFINED, Params(12, 100, 0.015e0)),
                        (CLASS_S, Params(12, 100, 0.015e0)),
                        (CLASS_W, Params(36, 400, 0.0015e0)),
                        (CLASS_A, Params(64, 400, 0.0015e0)),
                        (CLASS_B, Params(102, 400, 0.001e0)),
                        (CLASS_C, Params(162 ,400, 0.00067e0)),
                        (CLASS_D, Params(408, 500, 0.0003e0)),
                        (CLASS_E, Params(1020, 500, 0.0001e0)),
                        (CLASS_F, Params(2560, 500, 0.15e-4)),
                        ])

export CLASS