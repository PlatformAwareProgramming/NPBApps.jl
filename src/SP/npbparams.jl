# CLASS = A
#  
#  
#  This file is generated automatically by the setparams utility.
#  It sets the number of processors and the class of the NPB
#  in this directory. Do not modify it by hand.
#  

const MPIFC = "mpifc"
const FFLAGS = ""

const problem_size = 12
const niter_default = 100
const dt_default = 0.015
const convertdouble = false
const compiletime = "22 Nov 2023"
const npbversion="3.4.2"
const cs1="mpifort"
const cs2="$(MPIFC)"
const cs3="(none)" 
const cs4="(none)"
const cs5="-O3 -fallow-argument-mismatch -fopenmp"
const cs6="$(FFLAGS)"
const cs7="randi8"
