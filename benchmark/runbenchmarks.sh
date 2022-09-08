#!/bin/bash

# for benchmark in {"stencil_5th_order"}
# do
benchmark="stencil_5th_order"
    res_file=${benchmark}_results.csv
    rm -f $res_file
    echo "nthreads,t_flat,t_block\n" >> $res_file
    for nthreads in {1,2,4,6,8,12,16,24}
    do
        julia -t ${nthreads} -O3 --check-bounds=no ${benchmark}.jl
    done
# done