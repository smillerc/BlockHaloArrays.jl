#!/bin/bash

for benchmark in {"stencil_5th_order","stencil_2nd_order"}
do
    echo $benchmark.jl
    res_file=${benchmark}_results.csv
    rm -f $res_file
    echo "nthreads,t_flat_ns,t_block_ns" > $res_file
    for nthreads in {1,2,4,8,16}
    do
        julia --project=.. -t ${nthreads} -O3 --check-bounds=no ${benchmark}.jl
    done
done


