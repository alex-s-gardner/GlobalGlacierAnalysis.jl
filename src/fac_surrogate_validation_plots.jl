# Generate validation plots: Model vs Actual FAC as function of melt change
# Shows timeseries and scatter plots for visual inspection

using JLD2
using DecisionTree
using Statistics
using Random
using CairoMakie
using Printf

println("="^70)
println("GENERATING FAC SURROGATE VALIDATION PLOTS")
println("="^70)

# Load training data
cache_file = "/mnt/bylot-r3/data/gemb/raw/fac_surrogate_training_data_ratio.jld2"
println("\nLoading data from: $cache_file")
training_data = load(cache_file)

X_all = training_data["X"]
y_all = training_data["y"]  # Absolute FAC at scaled melt
feature_names = training_data["feature_names"]
lat_all = training_data["lat"]
lon_all = training_data["lon"]
elev_delta_all = training_data["elev_delta"]
pscale_all = training_data["pscale"]

println("  Loaded $(size(X_all, 1)) samples")

# Load model
model_file = "/mnt/bylot-r3/data/gemb/raw/fac_surrogate_random_forest.jld2"
println("Loading model from: $model_file")
fac_model_data = load(model_file)
model = fac_model_data["model"]

# Make predictions for all data
println("Generating predictions for all samples...")
y_pred_all = apply_forest(model, X_all)

# Extract features for analysis
melt_ratio_all = X_all[:, 1]  # melt_ratio
ref_fac_all = X_all[:, 6]     # ref_fac
temp_anomaly_all = X_all[:, 11]  # temp_anomaly

# Compute ΔFAC
Δfac_true = y_all .- ref_fac_all
Δfac_pred = y_pred_all .- ref_fac_all

println("\nCreating output directory...")
output_dir = "/mnt/bylot-r3/altim_figs/fac_surrogate"
mkpath(output_dir)

# ============================================================================
# PLOT 1: Model vs Actual FAC - Overall scatter
# ============================================================================
println("\nPlot 1: Overall model vs actual FAC scatter...")

fig1 = Figure(size=(1200, 1000))

# Panel A: Absolute FAC
ax1 = Axis(fig1[1, 1],
    xlabel="Actual FAC [m]",
    ylabel="Predicted FAC [m]",
    title="Model Performance: Absolute FAC",
    aspect=DataAspect())

# Density scatter with color
scatter!(ax1, y_all, y_pred_all,
    markersize=3, alpha=0.3, color=:blue)

# 1:1 line
xlims = (0, maximum(y_all))
lines!(ax1, [0, xlims[2]], [0, xlims[2]],
    color=:red, linestyle=:dash, linewidth=2, label="1:1 line")

# R² annotation
r2_abs = 1 - sum((y_all .- y_pred_all).^2) / sum((y_all .- mean(y_all)).^2)
rmse_abs = sqrt(mean((y_all .- y_pred_all).^2))
text!(ax1, 0.05, 0.95,
    text="R² = $(@sprintf("%.4f", r2_abs))\nRMSE = $(@sprintf("%.2f", rmse_abs)) m",
    space=:relative, align=(:left, :top), fontsize=14)

# Panel B: ΔFAC (Change)
ax2 = Axis(fig1[1, 2],
    xlabel="Actual ΔFAC [m]",
    ylabel="Predicted ΔFAC [m]",
    title="Model Performance: FAC Change (ΔFAC)",
    aspect=DataAspect())

scatter!(ax2, Δfac_true, Δfac_pred,
    markersize=3, alpha=0.3, color=:green)

# 1:1 line
xlims_delta = extrema(Δfac_true)
lines!(ax2, [xlims_delta[1], xlims_delta[2]], [xlims_delta[1], xlims_delta[2]],
    color=:red, linestyle=:dash, linewidth=2)

# R² annotation
r2_delta = 1 - sum((Δfac_true .- Δfac_pred).^2) / sum((Δfac_true .- mean(Δfac_true)).^2)
rmse_delta = sqrt(mean((Δfac_true .- Δfac_pred).^2))
text!(ax2, 0.05, 0.95,
    text="R² = $(@sprintf("%.4f", r2_delta))\nRMSE = $(@sprintf("%.2f", rmse_delta)) m",
    space=:relative, align=(:left, :top), fontsize=14)

# Panel C: Residuals vs melt_ratio
ax3 = Axis(fig1[2, 1],
    xlabel="Melt Ratio (mscale)",
    ylabel="Residual (Predicted - Actual) [m]",
    title="Residuals vs Melt Scaling")

residuals = Δfac_pred .- Δfac_true
scatter!(ax3, melt_ratio_all, residuals,
    markersize=3, alpha=0.2, color=:purple)
hlines!(ax3, [0], color=:red, linestyle=:dash, linewidth=2)

# Panel D: Residuals vs reference FAC
ax4 = Axis(fig1[2, 2],
    xlabel="Reference FAC [m]",
    ylabel="Residual (Predicted - Actual) [m]",
    title="Residuals vs Reference FAC")

scatter!(ax4, ref_fac_all, residuals,
    markersize=3, alpha=0.2, color=:orange)
hlines!(ax4, [0], color=:red, linestyle=:dash, linewidth=2)

save(joinpath(output_dir, "validation_overall_scatter.png"), fig1)
println("  Saved: validation_overall_scatter.png")

# ============================================================================
# PLOT 2: FAC vs Melt Ratio for Representative Locations
# ============================================================================
println("\nPlot 2: FAC vs melt ratio for representative locations...")

# Select 9 diverse locations (different latitude, elevation, climate)
unique_locs = unique([(lat_all[i], lon_all[i]) for i in 1:length(lat_all)])
Random.seed!(42)

# Try to get diverse sample
n_sample_locs = min(9, length(unique_locs))
sampled_locs = []

# Get locations at different latitudes
lat_bins = [0, 30, 60, 90]
for i in 1:length(lat_bins)-1
    locs_in_bin = [loc for loc in unique_locs if lat_bins[i] <= abs(loc[1]) < lat_bins[i+1]]
    if length(locs_in_bin) > 0
        push!(sampled_locs, rand(locs_in_bin))
        if length(sampled_locs) >= n_sample_locs
            break
        end
    end
end

# Fill remaining with random
while length(sampled_locs) < n_sample_locs && length(sampled_locs) < length(unique_locs)
    loc = rand(unique_locs)
    if !(loc in sampled_locs)
        push!(sampled_locs, loc)
    end
end

fig2 = Figure(size=(1600, 1200))

for (idx, loc) in enumerate(sampled_locs)
    row = div(idx - 1, 3) + 1
    col = mod(idx - 1, 3) + 1

    # Find all samples for this location
    loc_mask = [(lat_all[i], lon_all[i]) == loc for i in 1:length(lat_all)]

    if sum(loc_mask) == 0
        continue
    end

    # Extract data for this location
    melt_ratios = melt_ratio_all[loc_mask]
    fac_actual = y_all[loc_mask]
    fac_predicted = y_pred_all[loc_mask]
    ref_fac = ref_fac_all[loc_mask][1]  # Should be same for all at this location

    # Sort by melt_ratio for cleaner lines
    sort_idx = sortperm(melt_ratios)
    melt_ratios = melt_ratios[sort_idx]
    fac_actual = fac_actual[sort_idx]
    fac_predicted = fac_predicted[sort_idx]

    # Create subplot
    ax = Axis(fig2[row, col],
        xlabel="Melt Ratio (mscale)",
        ylabel="FAC [m]",
        title="Lat: $(@sprintf("%.1f", loc[1]))°, Lon: $(@sprintf("%.1f", loc[2]))°")

    # Reference line
    hlines!(ax, [ref_fac], color=:gray, linestyle=:dot, linewidth=2, label="Reference")

    # Actual values
    scatter!(ax, melt_ratios, fac_actual,
        markersize=8, color=:blue, label="Actual")

    # Predicted values
    scatter!(ax, melt_ratios, fac_predicted,
        markersize=8, color=:red, marker=:xcross, label="Predicted")

    # Add legend to first panel
    if idx == 1
        axislegend(ax, position=:rt, framevisible=true)
    end
end

save(joinpath(output_dir, "validation_fac_vs_melt_locations.png"), fig2)
println("  Saved: validation_fac_vs_melt_locations.png")

# ============================================================================
# PLOT 3: ΔFAC vs Melt Ratio binned
# ============================================================================
println("\nPlot 3: Binned ΔFAC vs melt ratio...")

# Bin by melt_ratio
melt_bins = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 6.0]
bin_centers = []
Δfac_true_means = []
Δfac_pred_means = []
Δfac_true_stds = []
Δfac_pred_stds = []

for i in 1:length(melt_bins)-1
    bin_mask = (melt_ratio_all .>= melt_bins[i]) .& (melt_ratio_all .< melt_bins[i+1])

    if sum(bin_mask) > 10  # Need reasonable sample size
        push!(bin_centers, (melt_bins[i] + melt_bins[i+1]) / 2)
        push!(Δfac_true_means, mean(Δfac_true[bin_mask]))
        push!(Δfac_pred_means, mean(Δfac_pred[bin_mask]))
        push!(Δfac_true_stds, std(Δfac_true[bin_mask]))
        push!(Δfac_pred_stds, std(Δfac_pred[bin_mask]))
    end
end

fig3 = Figure(size=(1200, 800))

ax1 = Axis(fig3[1, 1],
    xlabel="Melt Ratio (mscale)",
    ylabel="Mean ΔFAC [m]",
    title="FAC Change vs Melt Scaling (Binned)")

# Actual
errorbars!(ax1, bin_centers, Δfac_true_means, Δfac_true_stds,
    color=:blue, whiskerwidth=10)
scatter!(ax1, bin_centers, Δfac_true_means,
    markersize=12, color=:blue, label="Actual")

# Predicted
errorbars!(ax1, bin_centers, Δfac_pred_means, Δfac_pred_stds,
    color=:red, whiskerwidth=10)
scatter!(ax1, bin_centers, Δfac_pred_means,
    markersize=12, color=:red, marker=:xcross, label="Predicted")

hlines!(ax1, [0], color=:gray, linestyle=:dash)
vlines!(ax1, [1], color=:gray, linestyle=:dot, label="Reference (mscale=1)")

axislegend(ax1, position=:rb)

save(joinpath(output_dir, "validation_delta_fac_binned.png"), fig3)
println("  Saved: validation_delta_fac_binned.png")

# ============================================================================
# PLOT 4: Performance by Climate Zone
# ============================================================================
println("\nPlot 4: Performance by climate zone...")

# Define climate zones
lat_zones = [(0, 30, "Tropical"), (30, 60, "Mid-latitude"), (60, 90, "Polar")]
colors_zones = [:red, :orange, :blue]

fig4 = Figure(size=(1400, 1000))

# Panel A: Scatter by zone
ax1 = Axis(fig4[1, 1],
    xlabel="Actual ΔFAC [m]",
    ylabel="Predicted ΔFAC [m]",
    title="Model Performance by Climate Zone",
    aspect=DataAspect())

for (i, (lat_min, lat_max, zone_name)) in enumerate(lat_zones)
    zone_mask = (abs.(lat_all) .>= lat_min) .& (abs.(lat_all) .< lat_max)

    if sum(zone_mask) > 0
        scatter!(ax1, Δfac_true[zone_mask], Δfac_pred[zone_mask],
            markersize=4, alpha=0.3, color=colors_zones[i], label=zone_name)
    end
end

# 1:1 line
xlims_all = extrema(Δfac_true)
lines!(ax1, [xlims_all[1], xlims_all[2]], [xlims_all[1], xlims_all[2]],
    color=:black, linestyle=:dash, linewidth=2)

axislegend(ax1, position=:lt)

# Panel B: R² by zone
ax2 = Axis(fig4[1, 2],
    xlabel="Climate Zone",
    ylabel="R² on ΔFAC",
    title="R² by Climate Zone",
    xticks=(1:3, [z[3] for z in lat_zones]))

r2_by_zone = Float64[]
for (lat_min, lat_max, zone_name) in lat_zones
    zone_mask = (abs.(lat_all) .>= lat_min) .& (abs.(lat_all) .< lat_max)

    if sum(zone_mask) > 0
        Δfac_zone_true = Δfac_true[zone_mask]
        Δfac_zone_pred = Δfac_pred[zone_mask]
        r2_zone = 1 - sum((Δfac_zone_true .- Δfac_zone_pred).^2) /
                      sum((Δfac_zone_true .- mean(Δfac_zone_true)).^2)
        push!(r2_by_zone, r2_zone)
    else
        push!(r2_by_zone, NaN)
    end
end

barplot!(ax2, 1:3, r2_by_zone, color=colors_zones)
hlines!(ax2, [0.9], color=:red, linestyle=:dash, linewidth=2, label="Target (0.9)")
ylims!(ax2, 0, 1)

axislegend(ax2, position=:lb)

# Panel C: RMSE by zone
ax3 = Axis(fig4[2, 1],
    xlabel="Climate Zone",
    ylabel="RMSE on ΔFAC [m]",
    title="RMSE by Climate Zone",
    xticks=(1:3, [z[3] for z in lat_zones]))

rmse_by_zone = Float64[]
for (lat_min, lat_max, zone_name) in lat_zones
    zone_mask = (abs.(lat_all) .>= lat_min) .& (abs.(lat_all) .< lat_max)

    if sum(zone_mask) > 0
        Δfac_zone_true = Δfac_true[zone_mask]
        Δfac_zone_pred = Δfac_pred[zone_mask]
        rmse_zone = sqrt(mean((Δfac_zone_true .- Δfac_zone_pred).^2))
        push!(rmse_by_zone, rmse_zone)
    else
        push!(rmse_by_zone, NaN)
    end
end

barplot!(ax3, 1:3, rmse_by_zone, color=colors_zones)

# Panel D: Sample size by zone
ax4 = Axis(fig4[2, 2],
    xlabel="Climate Zone",
    ylabel="Number of Samples",
    title="Sample Distribution",
    xticks=(1:3, [z[3] for z in lat_zones]))

n_by_zone = Int[]
for (lat_min, lat_max, zone_name) in lat_zones
    zone_mask = (abs.(lat_all) .>= lat_min) .& (abs.(lat_all) .< lat_max)
    push!(n_by_zone, sum(zone_mask))
end

barplot!(ax4, 1:3, n_by_zone, color=colors_zones)

save(joinpath(output_dir, "validation_by_climate_zone.png"), fig4)
println("  Saved: validation_by_climate_zone.png")

# ============================================================================
# PLOT 5: Physical Constraint Validation
# ============================================================================
println("\nPlot 5: Physical constraint validation...")

fig5 = Figure(size=(1400, 1000))

# Panel A: FAC vs melt_ratio (should show negative trend)
ax1 = Axis(fig5[1, 1],
    xlabel="Melt Ratio (mscale)",
    ylabel="ΔFAC [m]",
    title="Physical Relationship: ΔFAC vs Melt Scaling")

# Bin and average to show trend
melt_fine_bins = range(minimum(melt_ratio_all), maximum(melt_ratio_all), length=50)
Δfac_trend_true = Float64[]
Δfac_trend_pred = Float64[]
melt_centers = Float64[]

for i in 1:length(melt_fine_bins)-1
    bin_mask = (melt_ratio_all .>= melt_fine_bins[i]) .& (melt_ratio_all .< melt_fine_bins[i+1])

    if sum(bin_mask) > 5
        push!(melt_centers, (melt_fine_bins[i] + melt_fine_bins[i+1]) / 2)
        push!(Δfac_trend_true, mean(Δfac_true[bin_mask]))
        push!(Δfac_trend_pred, mean(Δfac_pred[bin_mask]))
    end
end

lines!(ax1, melt_centers, Δfac_trend_true, color=:blue, linewidth=3, label="Actual")
lines!(ax1, melt_centers, Δfac_trend_pred, color=:red, linewidth=3, label="Predicted")
hlines!(ax1, [0], color=:gray, linestyle=:dash)
vlines!(ax1, [1], color=:gray, linestyle=:dot)

axislegend(ax1, position=:rb)

# Panel B: Check for negative predictions (should be minimal)
ax2 = Axis(fig5[1, 2],
    xlabel="Predicted FAC [m]",
    ylabel="Count",
    title="Distribution of Predicted FAC (check FAC ≥ 0)")

hist!(ax2, y_pred_all, bins=50, color=(:blue, 0.5))
vlines!(ax2, [0], color=:red, linewidth=3, label="FAC = 0")

n_negative = sum(y_pred_all .< 0)
text!(ax2, 0.95, 0.95,
    text="Negative predictions: $n_negative\n($(@sprintf("%.2f", 100*n_negative/length(y_pred_all)))%)",
    space=:relative, align=(:right, :top), fontsize=12)

# Panel C: Monotonicity check - FAC should decrease with melt at fixed location
ax3 = Axis(fig5[2, 1],
    xlabel="Rank by Melt Ratio",
    ylabel="ΔFAC [m]",
    title="Monotonicity Check (Sample Location)")

# Pick one location and show how ΔFAC changes with melt
if length(sampled_locs) > 0
    test_loc = sampled_locs[1]
    loc_mask = [(lat_all[i], lon_all[i]) == test_loc for i in 1:length(lat_all)]

    melt_loc = melt_ratio_all[loc_mask]
    Δfac_actual_loc = Δfac_true[loc_mask]
    Δfac_pred_loc = Δfac_pred[loc_mask]

    sort_idx = sortperm(melt_loc)

    lines!(ax3, 1:length(sort_idx), Δfac_actual_loc[sort_idx],
        color=:blue, linewidth=3, label="Actual", marker=:circle)
    lines!(ax3, 1:length(sort_idx), Δfac_pred_loc[sort_idx],
        color=:red, linewidth=3, label="Predicted", marker=:xcross)

    axislegend(ax3, position=:rb)
end

# Panel D: Prediction errors vs prediction magnitude
ax4 = Axis(fig5[2, 2],
    xlabel="|Predicted ΔFAC| [m]",
    ylabel="|Residual| [m]",
    title="Absolute Error vs Prediction Magnitude")

scatter!(ax4, abs.(Δfac_pred), abs.(residuals),
    markersize=3, alpha=0.2, color=:purple)

# Add trend line
abs_pred_bins = range(0, maximum(abs.(Δfac_pred)), length=20)
abs_error_trend = Float64[]
abs_pred_centers = Float64[]

for i in 1:length(abs_pred_bins)-1
    bin_mask = (abs.(Δfac_pred) .>= abs_pred_bins[i]) .& (abs.(Δfac_pred) .< abs_pred_bins[i+1])

    if sum(bin_mask) > 10
        push!(abs_pred_centers, (abs_pred_bins[i] + abs_pred_bins[i+1]) / 2)
        push!(abs_error_trend, mean(abs.(residuals[bin_mask])))
    end
end

lines!(ax4, abs_pred_centers, abs_error_trend, color=:red, linewidth=3)

save(joinpath(output_dir, "validation_physical_constraints.png"), fig5)
println("  Saved: validation_physical_constraints.png")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================
println("\n" * "="^70)
println("VALIDATION SUMMARY")
println("="^70)
println("\nOverall Performance:")
println("  R² (absolute FAC) = $(@sprintf("%.5f", r2_abs))")
println("  R² (ΔFAC) = $(@sprintf("%.5f", r2_delta)) ← PRIMARY METRIC")
println("  RMSE (ΔFAC) = $(@sprintf("%.3f", rmse_delta)) m")
println("  MAE (ΔFAC) = $(@sprintf("%.3f", mean(abs.(residuals)))) m")

println("\nPhysical Constraints:")
println("  Negative FAC predictions: $n_negative / $(length(y_pred_all)) ($(@sprintf("%.2f", 100*n_negative/length(y_pred_all)))%)")
println("  Mean ΔFAC at mscale=1: $(@sprintf("%.3f", mean(Δfac_true[melt_ratio_all .≈ 1.0]))) m (should be ~0)")

println("\nPerformance by Climate Zone:")
for (i, (lat_min, lat_max, zone_name)) in enumerate(lat_zones)
    println("  $zone_name (|lat| $lat_min-$lat_max°):")
    println("    R² = $(@sprintf("%.4f", r2_by_zone[i]))")
    println("    RMSE = $(@sprintf("%.3f", rmse_by_zone[i])) m")
    println("    N = $(n_by_zone[i]) samples")
end

println("\n" * "="^70)
println("PLOTS SAVED TO: $output_dir")
println("="^70)
println("  1. validation_overall_scatter.png - Overall performance scatter")
println("  2. validation_fac_vs_melt_locations.png - FAC vs melt for 9 locations")
println("  3. validation_delta_fac_binned.png - Binned ΔFAC vs melt ratio")
println("  4. validation_by_climate_zone.png - Performance by latitude zone")
println("  5. validation_physical_constraints.png - Physical validation checks")
println("\nReview these plots to visually inspect model accuracy!")
