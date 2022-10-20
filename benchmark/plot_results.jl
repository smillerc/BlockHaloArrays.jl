using CairoMakie
using CSV

function strong_scaling(P, N)
    S = 1 - P
    1 / (S + (P / N))
end

f = CSV.File("deluge_results.csv")
max_nthreads = maximum(f.nthreads)

fig = Figure()
ax = Axis(fig[1, 1], title="Strong Scaling",
    xlabel="Threads",
    ylabel="Speedup (vs 1 thread)",
)

lines!(ax, 1:max_nthreads, 1:max_nthreads, label="Linear Scaling", color=:black)
lines!(ax, 1:max_nthreads, strong_scaling.(0.99, 1:max_nthreads),
    label="99% Parallel", color=:green, linestyle=:dash)
lines!(ax, f.nthreads, first(f.t_block_ns) ./ f.t_block_ns, label="BlockHaloArray")
lines!(ax, f.nthreads, first(f.t_flat_ns) ./ f.t_flat_ns, label="@tturbo")

ax2 = Axis(fig[2, 1], title="Runtime vs Threads",
    xlabel="Threads", ylabel="Runtime[s]")

ax3 = Axis(fig[2, 1], yaxisposition=:right, ylabel="@tturbo / BlockHaloArray")
hidespines!(ax3)
hidexdecorations!(ax3)

lines!(ax2, f.nthreads, f.t_block_ns, label="BlockHaloArray")
lines!(ax2, f.nthreads, f.t_flat_ns, label="@tturbo")
lines!(ax3, f.nthreads, f.t_flat_ns ./ f.t_block_ns, color=:black)
lines!(ax2, [NaN], [NaN], label="speedup", color=:black)

fig[1, 2] = Legend(fig, ax,)
fig[2, 2] = Legend(fig, ax2,)

fig