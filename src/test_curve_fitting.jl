## define an model for which I would like to optimize the parameters to best fit data

# This function creates a "tilted" sine wave with offset and trend
eqn(t,p) =  p[1] .+ p[2] * t .+ p[3] * (1 / p[5]) * atan.(p[5] * sin.(2.0 * pi * (t .+ p[4])) ./ (1 .- p[5] * cos.(2.0 * pi * (t .+ p[4]))))

function eqn2(t,p)
    amp_scale = 1 / ((1 - (sqrt(1.0 - p[5]^2) + 2 * sqrt(1.0 - p[5]^4)) / 3) * (pi / 2 - 1) + 1)
    t_cycle = 2.0 * pi * (t .+ p[4])
    y = p[1] .+ p[2] * t .+ p[3] * (amp_scale / p[5]) * atan.(p[5] * sin.(t_cycle) ./ (1 .- p[5] * cos.(t_cycle)))
    return y
end
#    amp_scale = 1 / ((1 .- (sqrt.(1.0 .- p[5] .^ 2) .+ 2 * sqrt.(1.0 .- p[5] .^ 4)) / 3) * (pi / 2 - 1) .+ 1)


# time at which observations are made
t = collect(0:0.001:10) 

# model parameter values that I'd like to invert for given t and y
p = [0.0, 0.0, 1, 0.2, 0.9] 

# synthetic observations
y = eqn.(t, Ref(p)) .+ rand(length(t));
y2 = eqn2(t, p) .+ rand(length(t));
using Plots
scatter(t, y)

## Use optimization.jl
using Optimization, OptimizationPolyalgorithms, ForwardDiff
p0 = [0.1, 0.1, 1, 0.1, 0.1]

function cost(p0, _)
    sum(abs2, eqn.(t, Ref(p0)) .- y) / length(y)
end

using BenchmarkTools
@btime begin
    optf = OptimizationFunction(cost, Optimization.AutoForwardDiff())

    prob = OptimizationProblem(optf, p0)
    sol = solve(prob, PolyOpt())
end
plot!(t, eqn(t, sol.u))


# use LsqFit
using LsqFit
@btime fit = curve_fit($eqn, $t, $y, $p0)
fit = curve_fit(eqn, t, y, p0);
# 97.539 ms (2634 allocations: 98.18 MiB)
plot!(t, eqn(t, fit.param));

# add bounds for parameters
lb = [-20000., -2000., -5000., 0., 0.]
ub = [20000., 2000., 5000., 1., 1.]
@btime fit = curve_fit($eqn, $t, $y, $p0, lower=$lb, upper=$ub)
fit = curve_fit(eqn, t, y, p0, lower=lb, upper=ub);
#97.838 ms (2636 allocations: 98.18 MiB)
plot!(t, eqn(t, fit.param));

# 🌟🌟🌟🌟 add bounds for parameters & use :forwarddiff 🌟🌟🌟🌟
@btime fit = curve_fit($eqn2, $t, $y, $p0, lower=$lb, upper=$ub, autodiff=:forwarddiff)
fit = curve_fit(eqn2, t, y, p0, lower=lb, upper=ub);
# 42.311 ms (606 allocations: 54.59 MiB)
plot!(t, eqn(t, fit.param));

# add bounds for parameters & use :forwarddiff & normalize
lb = [-1., -5., -3., 0., 0.]
ub = [1., 5., 3., 1., 1.]
y0 = (y .- mean(y)) / std(y);
@btime fit = curve_fit($eqn, $t, $y0, $p0, lower=$lb, upper=$ub, autodiff=:forwarddiff)
fit = curve_fit(eqn, t, y, p0, lower=lb, upper=ub);
#44.509 ms (566 allocations: 64.20 MiB)
plot!(t, eqn(t, fit.param));

# use LsqFit
@btime fit = curve_fit($eqn, $t, $y, $p0;  autodiff=:forwarddiff)
fit = curve_fit(eqn, t, y, p0);
# 44.936 ms (564 allocations: 64.20 MiB)
plot!(t, eqn(t, fit.param))

# try eqn2
@btime fit = curve_fit($eqn2, $t, $y, $p0;  autodiff=:forwarddiff)
fit = curve_fit(eqn2, t, y, p0);
# 42.410 ms (603 allocations: 54.59 MiB)
plot!(t, eqn(t, fit.param))

# try eqn2
@btime fit = curve_fit($eqn2, $t, $y, $p0; autodiff=:finiteforward)
fit = curve_fit(eqn2, t, y, p0);
#  54.552 ms (1658 allocations: 47.54 MiB)
plot!(t, eqn(t, fit.param))



