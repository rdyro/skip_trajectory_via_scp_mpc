using PyPlot
using JLD2

function plot_single()
    @load "data_single2.jld2" Xactual Uactual levels N Xplan Uplan
    COLORS = ["#8c1515", "#2e2d29", "#175e54"]

    figure(figsize=(10,10))
    clf()
    plot(Xactual[1:3:end] .* 1e1, label=L"$v \: [m/s \: x \: 10^{-4}]$",
         color=COLORS[1], linewidth=3)
    plot(Xactual[2:3:end], label=L"$\theta \: [rad]$",
         color=COLORS[2], linewidth=3)
    plot(Xactual[3:3:end], label=L"$h \: [m \: x \: 10^{-5}]$",
         color=COLORS[3], linewidth=3)
    plot(Xplan[1:3:end] .* 1e1, label=L"$v_{plan} \: [m/s \: x \: 10^{-4}]$",
         color=COLORS[1], linewidth=3, linestyle="--")
    plot(Xplan[2:3:end], label=L"$\theta_{plan}\: [rad]$",
         color=COLORS[2], linewidth=3, linestyle="--")
    plot(Xplan[3:3:end], label=L"$h_{plan} \: [m \: x \: 10^{-5}]$",
         color=COLORS[3], linewidth=3, linestyle="--")
    legend(fontsize=15, loc=1)
    tight_layout()
    savefig("single-state.svg")

    figure(figsize=(10,10))
    plot(Uactual[1:end], label=L"$U$", color=COLORS[1], linewidth=3) 
    plot(Uplan[1:end], label=L"$U_{plan}$", color=COLORS[1], linewidth=3,
         linestyle="--")
    legend(fontsize=15, loc=1)
    tight_layout()
    savefig("single-control.svg")
    return


end

function plot_multi()
    @load "data_multi_not_noisey.jld2" Xactual Uactual levels N
    COLORS = ["#8c1515", "#2e2d29", "#175e54"]

    figure(figsize=(10,10))
    clf()
    plot(Xactual[1:3:end] .* 1e1, label=L"$v \: [m/s \: x \: 10^{-4}]$",
         color=COLORS[1], linewidth=3)
    plot(Xactual[2:3:end], label=L"$\theta \: [rad]$",
         color=COLORS[2], linewidth=3)
    plot(Xactual[3:3:end], label=L"$h \: [m \: x \: 10^{-5}]$",
         color=COLORS[3], linewidth=3)
    legend(fontsize=15, loc=1)
    tight_layout() 
    savefig("multi-state.svg")
    return

end

function plot_deviation()
    @load "deviation2.jld2" Xactual Uactual levels N Xplan Uplan
    COLORS = ["#8c1515", "#2e2d29", "#175e54"]

    figure(figsize=(10,6))
    clf()
    plot(Xactual[1:3:end] .* 1e1, label=L"$v \: [m/s \: x \: 10^{-4}]$",
         color=COLORS[1], linewidth=3)
    plot(Xactual[2:3:end], label=L"$\theta \: [rad]$",
         color=COLORS[2], linewidth=3)
    plot(Xactual[3:3:end], label=L"$h \: [m \: x \: 10^{-5}]$",
         color=COLORS[3], linewidth=3)
    plot(Xplan[1:3:end] .* 1e1, label=L"$v_{plan} \: [m/s \: x \: 10^{-4}]$",
         color=COLORS[1], linewidth=3, linestyle="--")
    plot(Xplan[2:3:end], label=L"$\theta_{plan}\: [rad]$",
         color=COLORS[2], linewidth=3, linestyle="--")
    plot(Xplan[3:3:end], label=L"$h_{plan} \: [m \: x \: 10^{-5}]$",
         color=COLORS[3], linewidth=3, linestyle="--")
    legend(fontsize=15, loc=1)
    tight_layout()
    savefig("deviation.svg")
    return

end

function plot_all()
    plot_single()
    plot_multi()
    plot_deviation()
    return
end

