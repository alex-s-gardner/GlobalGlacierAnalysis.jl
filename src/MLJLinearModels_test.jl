using MLJ
using MLJLinearModels
using Plots
using Random

# create data 
t = 1:0.01:10;
n = length(t);
gaussian_noise = randn(n) * 3;
outliers = rand((zeros(round(Int64, n/20))..., 6, -8, 100, -200, 178, -236, 77, -129, -50, -100, -45, -33, -114, -1929, -2000), n);

# measurment y
y = 10 .+ 10 * sin.(t) .+ 5 * t .+ gaussian_noise .+ outliers;

# design matrix 
X = hcat(ones(length(t)), sin.(t), t);

scale_penalty = false
fit_intercept = false

begin
    scatter(t, y; 
        markerstrokecolor=:match, 
        markerstrokewidth=0, 
        label = "observations", 
        ylim = (-70, 70),
        legend = :outerbottom,
        color = :grey,
        size = (700, 900)
    )

    # Base LSQ model fit
    println("Base Julia Linear Least Squares")
    @time θ = X \ y;
    plot!(t, X * θ, label="Base Julia Linear Least Squares", linewidth=2)

    regressor = LinearRegression(fit_intercept=fit_intercept);
    println(typeof(regressor))
    @time θ = fit(regressor, X, y)
    plot!(t, X * θ, label=typeof(regressor), linewidth = 2)

    regressor = HuberRegression(scale_penalty_with_samples=scale_penalty, fit_intercept=fit_intercept);
    println(typeof(regressor))
    @time θ = fit(regressor, X, y)
    plot!(t, X * θ, label=typeof(regressor), linewidth = 2)

    regressor = RidgeRegression(scale_penalty_with_samples=scale_penalty, fit_intercept=fit_intercept);
    println(typeof(regressor))
    @time θ = fit(regressor, X, y)
    plot!(t, X * θ, label=typeof(regressor), linewidth = 2)

    regressor = LassoRegression(scale_penalty_with_samples=scale_penalty, fit_intercept=fit_intercept);
    println(typeof(regressor))
    @time θ = fit(regressor, X, y)
    plot!(t, X * θ, label=typeof(regressor), linewidth = 2)

    regressor = ElasticNetRegression(scale_penalty_with_samples=scale_penalty, fit_intercept=fit_intercept);
    println(typeof(regressor))
    @time θ = fit(regressor, X, y)
    plot!(t, X * θ, label=typeof(regressor), linewidth = 2)

    regressor = QuantileRegression(scale_penalty_with_samples=scale_penalty, fit_intercept=fit_intercept);
    println(typeof(regressor))
    @time θ = fit(regressor, X, y)
    plot!(t, X * θ, label=typeof(regressor), linewidth = 2)

    regressor = LADRegression(scale_penalty_with_samples=scale_penalty, fit_intercept=fit_intercept);
    println(typeof(regressor))
    @time θ = fit(regressor, X, y)
    plot!(t, X * θ, label=typeof(regressor), linewidth = 2)

    regressor = GeneralizedLinearRegression(scale_penalty_with_samples=scale_penalty, fit_intercept=fit_intercept);
    println(typeof(regressor))
    @time θ = fit(regressor, X, y)
    plot!(t, X * θ, label=typeof(regressor), linewidth = 2)

    regressor = RobustRegression(scale_penalty_with_samples=scale_penalty, fit_intercept=fit_intercept);
    println(typeof(regressor))
    @time θ = fit(regressor, X, y)
    plot!(t, X * θ, label=typeof(regressor), linewidth = 2)
end