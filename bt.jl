using Pkg
Pkg.activate(".")
using NPBApps
#BT.go(BT.CLASS_W)
BT.go(Array{Int64}([12,12,6]), 60, 0.010)
