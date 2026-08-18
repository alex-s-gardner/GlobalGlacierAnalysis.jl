# Validate R² on FAC CHANGE (ΔFAC), not absolute FAC
# Critical: Model should predict the perturbation response, not just absolute values

using JLD2
using DecisionTree
using Statistics
using Random

println("="^70)
println("VALIDATING R² ON FAC CHANGE (ΔFAC)")
println("="^70)

# Load training data
cache_file = "/mnt/bylot-r3/data/gemb/raw/fac_surrogate_training_data_ratio.jld2"
println("\nLoading data from: $cache_file")
training_data = load(cache_file)

X_all = training_data["X"]
y_all = training_data["y"]  # This is absolute FAC at scaled melt
feature_names = training_data["feature_names"]
lat_all = training_data["lat"]
lon_all = training_data["lon"]

println("  Loaded $(size(X_all, 1)) samples")

# Subset to 1000 locations (same as training)
use_subset = true
n_locations_subset = 1000

if use_subset
    loc_keys = [(lat_all[i], lon_all[i]) for i in 1:length(lat_all)]
    unique_locs = unique(loc_keys)
    Random.seed!(42)
    selected_locs = Set(shuffle(unique_locs)[1:min(n_locations_subset, length(unique_locs))])
    subset_mask = [loc_keys[i] in selected_locs for i in 1:length(lat_all)]

    X_all = X_all[subset_mask, :]
    y_all = y_all[subset_mask]
    lat_all = lat_all[subset_mask]
    lon_all = lon_all[subset_mask]
end

# Location-stratified split
loc_keys = [(lat_all[i], lon_all[i]) for i in 1:length(lat_all)]
unique_locs = unique(loc_keys)
n_locs = length(unique_locs)

Random.seed!(42)
shuffled_locs = shuffle(unique_locs)
n_train_locs = round(Int, 0.70 * n_locs)
n_val_locs = round(Int, 0.15 * n_locs)

train_locs = Set(shuffled_locs[1:n_train_locs])
val_locs = Set(shuffled_locs[n_train_locs+1:n_train_locs+n_val_locs])
test_locs = Set(shuffled_locs[n_train_locs+n_val_locs+1:end])

train_mask = [loc_keys[i] in train_locs for i in 1:length(lat_all)]
test_mask = [loc_keys[i] in test_locs for i in 1:length(lat_all)]

X_train = X_all[train_mask, :]
y_train = y_all[train_mask]
X_test = X_all[test_mask, :]
y_test = y_all[test_mask]

println("\nDataset splits:")
println("  Train: $(size(X_train, 1)) samples")
println("  Test: $(size(X_test, 1)) samples")

# Load model
model_file = "/mnt/bylot-r3/data/gemb/raw/fac_surrogate_random_forest.jld2"
println("\nLoading model from: $model_file")
fac_model_data = load(model_file)
model = fac_model_data["model"]

# Make predictions
println("\nGenerating predictions...")
y_train_pred = apply_forest(model, X_train)
y_test_pred = apply_forest(model, X_test)

# Extract reference FAC from features (column 6: ref_fac)
ref_fac_train = X_train[:, 6]
ref_fac_test = X_test[:, 6]

# Compute CHANGE (ΔFAC = FAC_scaled - FAC_ref)
Δfac_train_true = y_train .- ref_fac_train
Δfac_train_pred = y_train_pred .- ref_fac_train

Δfac_test_true = y_test .- ref_fac_test
Δfac_test_pred = y_test_pred .- ref_fac_test

# R² on ABSOLUTE FAC (what I reported before)
function r2_score(y_true, y_pred)
    ss_res = sum((y_true .- y_pred).^2)
    ss_tot = sum((y_true .- mean(y_true)).^2)
    return 1 - ss_res / ss_tot
end

r2_train_absolute = r2_score(y_train, y_train_pred)
r2_test_absolute = r2_score(y_test, y_test_pred)

# R² on CHANGE (ΔFAC) - THE CRITICAL METRIC
r2_train_change = r2_score(Δfac_train_true, Δfac_train_pred)
r2_test_change = r2_score(Δfac_test_true, Δfac_test_pred)

# RMSE on change
rmse_train_change = sqrt(mean((Δfac_train_true .- Δfac_train_pred).^2))
rmse_test_change = sqrt(mean((Δfac_test_true .- Δfac_test_pred).^2))

# MAE on change
mae_train_change = mean(abs.(Δfac_train_true .- Δfac_train_pred))
mae_test_change = mean(abs.(Δfac_test_true .- Δfac_test_pred))

println("\n" * "="^70)
println("RESULTS: R² ON ABSOLUTE FAC (what I reported)")
println("="^70)
println("  Train R² = $(round(r2_train_absolute, digits=5))")
println("  Test R² = $(round(r2_test_absolute, digits=5))")

println("\n" * "="^70)
println("RESULTS: R² ON CHANGE (ΔFAC) - THE REAL TEST")
println("="^70)
println("  Train R² = $(round(r2_train_change, digits=5))")
println("  Test R² = $(round(r2_test_change, digits=5))")
println("  Test RMSE = $(round(rmse_test_change, digits=3)) m")
println("  Test MAE = $(round(mae_test_change, digits=3)) m")

println("\n" * "="^70)
println("COMPARISON")
println("="^70)
println("Metric              | Absolute FAC | ΔFAC (Change)")
println("--------------------|--------------|---------------")
println("Test R²             | $(rpad(round(r2_test_absolute, digits=5), 12)) | $(round(r2_test_change, digits=5))")
println("Test RMSE           | $(rpad(round(sqrt(mean((y_test .- y_test_pred).^2)), digits=3), 12)) m | $(round(rmse_test_change, digits=3)) m")
println("Test MAE            | $(rpad(round(mean(abs.(y_test .- y_test_pred)), digits=3), 12)) m | $(round(mae_test_change, digits=3)) m")

println("\n" * "="^70)
println("INTERPRETATION")
println("="^70)

if r2_test_change >= 0.9
    println("✓ SUCCESS: R² on ΔFAC = $(round(r2_test_change, digits=5)) ≥ 0.9")
    println("  Model accurately predicts FAC change due to melt scaling")
else
    println("✗ FAILURE: R² on ΔFAC = $(round(r2_test_change, digits=5)) < 0.9")
    println("  Model may be mostly predicting ref_fac, not the actual change")
    println("  Need to retrain with ΔFAC as target, not absolute FAC")
end

println("\nΔFAC statistics (test set):")
println("  Mean: $(round(mean(Δfac_test_true), digits=3)) m")
println("  Std: $(round(std(Δfac_test_true), digits=3)) m")
println("  Range: [$(round(minimum(Δfac_test_true), digits=2)), $(round(maximum(Δfac_test_true), digits=2))] m")

println("\n" * "="^70)
println("RECOMMENDATION")
println("="^70)

if r2_test_change < 0.9
    println("Need to retrain model with ΔFAC as target:")
    println("  y_all = fac_scaled_vec - ref_fac_vec  # Train on change, not absolute")
    println("  This forces model to learn perturbation response, not memorize ref_fac")
else
    println("Current model is adequate - accurately predicts FAC changes")
end
