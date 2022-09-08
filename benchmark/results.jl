using CairoMakie
using CSV


f = CSV.File("results_5th_order.csv")
max_nthreads = maximum(f.nthreads)

fig = Figure()
ax = Axis(fig[1, 1], title="Speedup",
    xlabel="Threads",
    ylabel="Speedup (vs 1 thread)",
)


lines!(ax, 1:max_nthreads, 1:max_nthreads, label="Linear Scaling", color=:black)
lines!(ax, f.nthreads, first(f.t_block) ./ f.t_block, label="BlockHaloArray")
lines!(ax, f.nthreads, first(f.t_flat) ./ f.t_flat, label="@tturbo")

ax2 = Axis(fig[2, 1], title="Strong Scaling",
    xlabel="Threads",
    # xticks=1:10,
    ylabel="Runtime[s]")

lines!(ax2, f.nthreads, f.t_block, label="BlockHaloArray")
lines!(ax2, f.nthreads, f.t_flat, label="@tturbo")

fig[1, 2] = Legend(fig, ax,)

fig