using Pkg
Pkg.activate(".")
using NPBApps
#LU.go(LU.CLASS_S)
LU.go(12, 12, 12, 50, 50, 0.5)
