using .Threads, LoopVectorization, ThreadPools, ThreadPinning

function _blocked(A, blockid)
    ilo, ihi, jlo, jhi = A.loop_limits[blockid]
    data = @views A.blocks[blockid]
    @turbo for j in jlo:jhi
        for i in ilo:ihi
            data[i, j] = i * j + 0.25(
                data[i-2, j] + data[i-1, j] + data[i+2, j] + data[i+1, j] +
                data[i, j] +
                data[i, j-2] + data[i, j-1] + data[i, j+2] + data[i, j+1])
        end
    end
end

function blocked_version(A)
    @sync for tid in 1:nthreads()
        ThreadPools.@tspawnat tid _blocked(A, tid)
    end
end

function tturbo_version(data, nhalo)
    ilo = first(axes(data, 1)) + nhalo
    ihi = last(axes(data, 1)) - nhalo
    jlo = first(axes(data, 2)) + nhalo
    jhi = last(axes(data, 2)) - nhalo

    @tturbo for j in jlo:jhi
        for i in ilo:ihi
            data[i, j] = i * j + 0.25(
                data[i-2, j] + data[i-1, j] + data[i+2, j] + data[i+1, j] +
                data[i, j] +
                data[i, j-2] + data[i, j-1] + data[i, j+2] + data[i, j+1])
        end
    end
end

pinthreads(:compact)

ni = 5000
nj = 5000
nhalo = 4
A_block = BlockHaloArray((ni, nj), nhalo);
A_flat = rand(ni + 2nhalo, nj + 2nhalo);

tturbo_version(A_flat, nhalo)
blocked_version(A_block)

t_elaps = 0
max_iter = 1000
for iter in 1:max_iter
    global t_elaps += @elapsed tturbo_version(A_flat, nhalo)
end
t_tturbo = t_elaps / max_iter

t_elaps = 0
for iter in 1:max_iter
    global t_elaps += @elapsed blocked_version(A_block)
end
t_blocked = t_elaps / max_iter

@show t_blocked t_tturbo
@show t_tturbo/ t_blocked

nothing

