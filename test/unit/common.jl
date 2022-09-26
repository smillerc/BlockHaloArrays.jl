using .Threads
using BlockHaloArrays
using EllipsisNotation
using ThreadPools: @tspawnat
using Test
# using TestItems, TestItemRunner

function init(A, thread_id)
    dom = domainview(A, thread_id)
    fill!(dom, thread_id)
end
