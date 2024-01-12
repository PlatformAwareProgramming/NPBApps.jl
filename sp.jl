using Pkg
Pkg.activate(".")
using NPBApps
#SP.go(SP.CLASS_C)
SP.go(Array{Int64}([12,12,12]), 100, 0.015)
