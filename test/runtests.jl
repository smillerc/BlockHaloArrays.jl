using TestItemRunner

@run_package_tests




# @testitem "Flat index testing" begin

#     dims = (4, 50, 50)
#     global_idx = (1, 26, 26)
#     halodims = (2, 3)
#     A = BlockHaloArray(dims, halodims, 2, 2)

#     block_idx = ntuple(i ->
#             BlockHaloArrays.get_block_idx(global_idx[A.halodims[i]],
#                 A._cummulative_blocksize_per_dim,
#                 A.halodims[i],
#                 A.halodims),
#         length(A.halodims))

#     local_idx = ntuple(i ->
#             BlockHaloArrays.get_local_idx(global_idx,
#                 block_idx,
#                 A.nhalo,
#                 A._cummulative_blocksize_per_dim,
#                 i,
#                 A.halodims),
#         length(A.globaldims))

#     @test block_idx == (2, 1)
#     @test local_idx == (1, 3, 28)
# end
