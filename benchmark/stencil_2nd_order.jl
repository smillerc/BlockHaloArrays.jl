import Pkg

Pkg.activate("..")

include("stencils.jl")

recon = MUSCLReconstruction(MinMod())
nhalo = 2
run_stencil_benchmark(recon, nhalo, Float64, "stencil_2nd_order_results.csv")
