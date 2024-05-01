using Altim
using Plots
using Extents
using Statistics
using Flux
using OptimizationOptimJL
using Extents
using BenchmarkTools
using BinStatistics
using FastGeoProjections

geotile_width = 2; #geotile width [degrees]
grid = (node_spacing = 100, node_width = 500); # [(node_spacing = 500, node_width = 500 * 2)]

project_id = :v01;
domain = :landice;
#domain = :all;

paths = project_paths(project_id = project_id);
products = project_products(; project_id=:v01);
geotiles = project_geotiles(; geotile_width = geotile_width, domain = domain);


dems = [:cop30_v2]; 
# missions = [:gedi, :icesat, :icesat2, :hugonnet]
# --------------------------------------------------------------------------
ext = Extent(X=(-126.9, -126.1), Y=(51.1, 51.8))
geotiles = geotile_subset!(geotiles, ext);

geotile = first(geotiles)
dem = first(dems)

@time df = geotile_merge_height(geotile, dem, paths; products=products)


dmad = df.mad_ref .- df.mad_offset;
valid = .!(df.dx .== 0);
histogram(dmad[valid], bins = -5:.1:5)
#geotilefn = joinpath(geotile_dir, geotile[:id] * ".arrow");
 

# create a regularly spaced grid using a local low distortion projection
df, epsg = geotile_utm!(df)

dh0 = df.height .- df.height_reference;

## USING MLJ
using MLJ
using MLJLinearModels
scale_penalty = false
fit_intercept = false

function elevation_change_linear_model(t, z; t_intercept=2015)
    X = hcat(ones(size(t)), (t .- t_intercept), (t .- t_intercept).^2, cos.(2 * pi * t), sin.(2 * pi * t),  (z .- mean(z)))
    return X
end

t0 = 2000:0.1:2023;
dz0 = zeros(size(t0));
X0 = elevation_change_linear_model(t0, dz0);


d = 100;
xpt = 693023;
ypt = 5707318;

d = 100
xpt = 699899
ypt = 5704911

d = 100
xpt = 692964
ypt = 5708435

xy = Tuple{Float64,Float64}[]
for i in eachindex(df.X)
    push!(xy,(df.X[i], df.Y[i]))
end

d = 250
begin
    (xpt, ypt) = rand(xy, 1)[1]

    ext = Extent(X = (xpt - d, xpt + d), Y = (ypt - d, ypt + d))


    ext = Extent{(:X, :Y),Tuple{Tuple{Float32,Float32},Tuple{Float32,Float32}}}((X=(-1.4039834f7, -1.4038461f7), Y=(6.7022205f6, 6.702905f6)))
    trans = FastGeoProjections.Transformation("EPSG:3857", epsg; always_xy=true);
    XY1 = trans.(ext.X,ext.Y)
    ext = Extent(X=(min(XY1[1][1][1], XY1[2][1][1]), max(XY1[1][1][1], XY1[2][1][1])), Y=(min(XY1[1][2][1], XY1[2][2][1]), max(XY1[1][2][1], XY1[2][2][1])))
   
    ind = within.(Ref(ext), df.X, df.Y)  .& .!isnan.(dh0) 

    t = df.decyear[ind]
    z = df.height_reference[ind];
    y = dh0[ind]
    w = 1 ./ df.error[ind];
    w = ones(size(w))
    product = df.product[ind]

    X = elevation_change_linear_model(t,z);

    up = unique(product)
    f = product .== up[1]
    scatter(t[f], y[f]; markerstrokecolor=:match, markerstrokewidth=0, label=up[1], ylims = (-1000, 1000))
    if length(up) >1
        for u in up[2:end]
            f = product .== u
                scatter!(t[f], y[f]; markerstrokecolor=:match, markerstrokewidth=0, label=u)
        end
    end
    
  
    # Base LSQ model fit
    println("Base LSQ")
    @time θ = (X.*w) \ (y.*w);
    plot!(t0, X0 * θ, label="Base lsq")

    
    # index X[:,2:end] as fit includes offset by default
    println("LinearRegression")
    θ = fit(LinearRegression(fit_intercept=fit_intercept), X.*w, y.*w)
    plot!(t0, X0 * θ,label="linear")

    # Huber produces results similar to Robust but is slighlty faster. 
    println("HuberRegression")
    θ = fit(HuberRegression(fit_intercept=fit_intercept, scale_penalty_with_samples = scale_penalty), X.*w, y.*w)
    plot!(t0, X0 * θ, label="huber", linestyle=:dash)

    println("RobustRegression")
    θ = fit(RobustRegression(fit_intercept=fit_intercept, scale_penalty_with_samples=scale_penalty), X.*w, y.*w)
    plot!(t0, X0 * θ, label="robust", linestyle=:dash)

    #=
    println("RidgeRegression")
    @time θ = fit(RidgeRegression(fit_intercept=fit_intercept, scale_penalty_with_samples=scale_penalty), X.*w, y.*w)
    plot!(t0, X0 * θ, label="ridge")

    println("LassoRegression")
    @time θ = fit(LassoRegression(fit_intercept=fit_intercept, scale_penalty_with_samples=scale_penalty), X.*w, y.*w)
    plot!(t0, X0 * θ, label="lasso")

    println("ElasticNetRegression")
    @time θ = fit(ElasticNetRegression(fit_intercept=fit_intercept, scale_penalty_with_samples=scale_penalty), X.*w, y.*w)
    plot!(t0, X0 * θ, label="elastic")

    println("QuantileRegression")
    @time θ = fit(QuantileRegression(fit_intercept=fit_intercept, scale_penalty_with_samples=scale_penalty), X.*w, y.*w)
    plot!(t0, X0 * θ, label="quantile")

    println("LADRegression")
    @time θ = fit(LADRegression(fit_intercept=fit_intercept, scale_penalty_with_samples=scale_penalty), X.*w, y.*w)
    plot!(t0, X0 * θ, label="lad")

    println("GeneralizedLinearRegression")
    @time θ = fit(GeneralizedLinearRegression(fit_intercept=fit_intercept, scale_penalty_with_samples=scale_penalty), X.*w, y.*w)
    plot!(t0, X0 * θ, label="generalized")
    =#

end