#!/bin/bash

# for nthreads in {2..4..2}
# do
#     julia -t ${nthreads} -O3 --check-bounds=no stencil.jl
# done

for nthreads in {6..8..2}
do
    julia -t ${nthreads} -O3 --check-bounds=no stencil_5th_order.jl
done
