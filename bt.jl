using Pkg
Pkg.activate(".")
using NPBApps
#BT.go(BT.CLASS_S)
BT.go(Array{Int64}([12,12,12]), 60, 0.010)