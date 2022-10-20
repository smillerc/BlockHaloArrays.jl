# BlockHaloArrays

The `BlockHaloArray` type is an array-like type (does not extend `AbstractArray`) that is designed for shared-memory multi-threaded workloads. Typical shared-memory parallelization is done via multi-threaded loops (with `@threads`, `@tturbo`, etc.), i.e.:
```julia
Threads.@threads for j in jlo:jhi
    for i in ilo:ihi
        A[i,j] = ...
    end
end
```
However, this type of parallelization tends not to scale well, especially for memory-bandwidth limited workloads. This is especially true for stencil operations, where functions access memory at addresses such as `A[i,j], A[i+1,j], A[i-1,j+1]`. The aim of the `BlockHaloArray` type do mini-domain decomposition into separate thread-specific chunks, or blocks, of memory for each thread to operate on without data race conditions or memory bandwidth contention. Because the context is shared memory however, certain operations are more efficient compared to the typical MPI domain decomposition method.Each thread will be able to operate on it's own block of memory and maintain NUMA-aware locality and communication is done via copies, rather than MPI library calls. Future releases will include an MPI-aware `BlockHaloArray` that will facilitate the MPI+Threads hybrid parallelization.

An example of a blocked-style stencil operation looks like the following. Note that it can still take advantage of the single-threaded `@turbo` macro to vectorize the nested-loop.
```julia
function stencil!(A::BlockHaloArray, blockid)
    ilo, ihi, jlo, jhi = A.loop_limits[blockid]
    data = @views A.blocks[blockid]
    @turbo for j in jlo:jhi
        for i in ilo:ihi
            data[i, j] = i * j + 0.25(data[i-2, j] + data[i-1, j] + data[i+2, j] + data[i+1, j] +
                                      data[i, j] +
                                      data[i, j-2] + data[i, j-1] + data[i, j+2] + data[i, j+1])
        end
    end
end
```

Then the `stencil!()` function is split among threads via:
```julia
@sync for block_id in 1:nthreads()
    @tspawnat block_id stencil!(A, block_id)
end
```

`BlockHaloArray` types should be used in a task-based parallel workload, which has better scaling than multi-threaded loops. The only synchronization required is during the exchange of halo data.

Benchmark results are coming soon...

## Exports

### Types

- `BlockHaloArray`: A blocked array-like type (does not extent `AbstractArray`) to be used in a shared-memory type system. This partitions an `Array` into N blocks that are to be operated on by threads.

#### Constructors
```julia
BlockHaloArray(dims::NTuple{N,Int}, halodims::NTuple{N2,Int}, nhalo::Integer, nblocks=nthreads(); T=Float64, use_numa=true)
BlockHaloArray(dims::NTuple{N,Int}, nhalo::Integer, nblocks=nthreads(); T=Float64, use_numa=true)
BlockHaloArray(A::AbstractArray{T,N}, nhalo::Integer, nblocks=nthreads())
BlockHaloArray(A::AbstractArray{T,N}, halodims::NTuple{N2,Integer}, nhalo::Integer, nblocks=nthreads()) 
```

#### Example
```julia
using .Threads
using ThreadPools
using BlockHaloArrays

dims = (30, 20) # a 30 x 20 matrix
nhalo = 2 # number of halo entries in each dimension
nblocks = nthreads() # nblocks must be ≤ nthreads() or a warning will be issued

A = BlockHaloArray(dims, nhalo, nblocks; T=Float64) # create with empty data
B = BlockHaloArray(rand(dims...), nhalo, nblocks) # create from an existing array

# Create a 3D array, but only do halo exchange on the 2nd and 3rd dimensions
newdims = (4, 100, 200)
haloax = (2, 3) # these must be monotonically increasing and the outer-most dims, i.e. can't be (1, 2), or (1, 3)
C = BlockHaloArray(newdims, haloax, nhalo, nblocks; T=Float64)

# Fill each block with it's corresponding threadid value
function init(A)
    dom = domainview(A, threadid())
    fill!(dom, threadid())
end

# Let each thread initialize it's own block
@sync for tid in 1:nblocks(A)
    @tspawnat tid init(A1)
end

# Let each block sync its halo region with its neighbors
updatehalo!(A)

Anew = flatten(A) # Anew is a 2D Array -- no longer a BlockHaloArray
```

- `MPIBlockHaloArray`: (to be implemented) An MPI-aware version of a `BlockHaloArray`. This facilitates halo communication between nodes via MPI.

### Functions

 - `flatten(A)`: Return a flattened `Array` of `A`. This is a copy.
 - `repartition(A, nblocks)`: Repartition the BlockHaloArray into a different block layout 
 - `updatehalo!(A, include_periodic_bc=false)`: Sync all block halo regions. This uses a `@sync` loop with @tspawnat to sync each block.
 - `updateblockhalo!(A, block_id, include_periodic_bc=false)`: Sync the halo regions of a single block. No `@sync` or `@spawn` calls are used.
 - `domainview(A, blockid)`: Return a view of a the domain (no halo regions) of block `blockid`
 - `onboundary(A, blockid)`: Is the current block `blockid` on a boundary? This is used help apply boundary conditions in a user code.

Additional overloaded methods include
```julia
eltype(A::AbstractBlockHaloArray) = eltype(first(A.blocks))

size(A::AbstractBlockHaloArray) = A.globaldims
size(A::AbstractBlockHaloArray, dim) = size.(A.blocks, dim)

axes(A::AbstractBlockHaloArray) = axes.(A.blocks)
axes(A::AbstractBlockHaloArray, dim) = axes.(A.blocks, dim)

first(A::AbstractBlockHaloArray) = first(A.blocks)
firstindex(A::AbstractBlockHaloArray) = firstindex(A.blocks)

last(A::AbstractBlockHaloArray) = last(A.blocks)
lastindex(A::AbstractBlockHaloArray) = lastindex(A.blocks)

eachindex(A::AbstractBlockHaloArray) = eachindex(A.blocks)

nblocks(A::AbstractBlockHaloArray) = length(A.blocks)
```
