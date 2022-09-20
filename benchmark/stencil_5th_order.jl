import Pkg

Pkg.activate("..")

include("stencils.jl")

recon = MLP5Reconstruction()
nhalo = 3
run_stencil_benchmark(recon, nhalo, Float64, "stencil_5th_order_results.csv")