using .Threads
using BlockHaloArrays
using EllipsisNotation
using BlockHaloArrays: @tspawnat
using Test
# using JET
# using TestItems, TestItemRunner

function init(A, thread_id)
    dom = domainview(A, thread_id)
    fill!(dom, thread_id)
end
