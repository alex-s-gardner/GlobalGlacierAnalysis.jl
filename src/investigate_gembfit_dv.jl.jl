 



using DataFrames
using GeoDataFrames
geotile_width =2;

single_geotile_test = "lat[+58+60]lon[-138-136]"; #GGA.geotiles_golden_test[7]

path2runs = path2runs_synthesized

begin
fig = Figure();
ax = Axis(fig[1, 1]; xlabel="date", ylabel="dv [km^3]", title="$single_geotile_test");
#for file_2_select in 1:length(path2runs[1:5:end])
file_2_select = 1
binned_synthesized_file = path2runs[file_2_select];

    params = GGA.binned_filled_fileparts(binned_synthesized_file);

    file_date = Dates.unix2datetime(mtime(binned_synthesized_file));

    dh = FileIO.load(binned_synthesized_file, "dh_hyps");

    dgeotile = dims(dh, :geotile);
    area_km2 = GGA._geotile_area_km2(; params.surface_mask, geotile_width)[geotile=At(val(dgeotile))];

    # Convert elevation change to volume change
    dv_altim = GGA.dh2dv_geotile(dh,area_km2);

    
    lines!(ax, dv_altim[geotile = At(single_geotile_test)]);

    idx = .!isnan.(dv_altim[geotile = At(single_geotile_test)]);
    println("mean = $(round(mean(dv_altim[geotile = At(single_geotile_test)][idx]), digits=2)), file = $binned_synthesized_file, date = $file_date");


display(fig)



synthesized_gemb_fit = replace(binned_synthesized_file, ".jld2" => "_gembfit.arrow");
gemb_fit = GeoDataFrames.read(synthesized_gemb_fit);

fit_file_date = Dates.unix2datetime(mtime(synthesized_gemb_fit))

gt_index = findfirst(geotiles.id .== single_geotile_test);
println("raw - pscale = $(gemb_fit[gt_index, :pscale]), ΔT = $(gemb_fit[gt_index, :ΔT])");

gemb_fit[gt_index, :pscale];
gemb_fit[gt_index, :ΔT];
binned_synthesized_dv_file = replace(binned_synthesized_file, ".jld2" => "_gembfit_dv.jld2");
geotiles = GGA.FileIO.load(binned_synthesized_dv_file, "geotiles");

# group by rgi and sum
geotiles_reg = groupby(geotiles, :rgiid);

# sum across timeseries
vars2sum = ["dv", "dm", "dv_altim", "dm_altim"];
regions0 = DataFrames.combine(geotiles_reg, vars2sum .=> Ref ∘ sum; renamecols=false);

lines(regions0[findfirst(regions0.rgiid .== 1),:dv]);


geotile_row = geotiles[findfirst(geotiles.id .== single_geotile_test), :];

#binned_synthesized_file  = path2runs[1]
#binned_synthesized_dv_file = replace(binned_synthesized_file, ".jld2" => "_gembfit_dv.jld2")
#geotiles = GGA.FileIO.load(binned_synthesized_dv_file, "geotiles")
#gt_index = findfirst(geotiles.id .== single_geotile_test)
#geotile_row = geotiles[gt_index, :];
gemb_date = colmetadata(geotile_row, "dm")["date"];
altim_date = colmetadata(geotile_row, "dm_altim")["date"];


fig = Figure();;
ax = Axis(fig[1, 1]);
lines!(ax, gemb_date, geotile_row.dv; label="dv");
lines!(ax, altim_date, geotile_row.dv_altim; label="dv_altim");
axislegend(ax, position=:rt);
display(fig);

geotile_row.pscale
geotile_row.ΔT

println("final - pscale = $(geotile_row.pscale), ΔT = $(geotile_row.ΔT)");

end

gt_index = findfirst(geotiles.id .== single_geotile_test)

geotile_row = geotiles[gt_index, :];

pscale = round(geotile_row["pscale"], digits=2)
ΔT = round(geotile_row["ΔT"], digits=2)

println("pscale = $pscale, ΔT = $ΔT, file = $(binned_synthesized_dv_file), last modified: $(Dates.unix2datetime(mtime(binned_synthesized_dv_file)))")





gemb_index = (gemb_date .> DateTime(2010, 1, 1)) .& (gemb_date .< DateTime(2015, 1, 1));
altim_index = (altim_date .> DateTime(2010, 1, 1)) .& (altim_date .< DateTime(2015, 1, 1));

#pscale = geotile_row["pscale"];
#ΔT = geotile_row["ΔT"];

begin
    varname = "dv";
    f = Figure();
    ax = Axis(f[1, 1], title="$single_geotile_test: $varname [pscale = $pscale, ΔT = $ΔT]");
    lines!(ax, gemb_date, geotile_row[varname] .- mean(geotile_row[varname][gemb_index]), label="GEMB");
    lines!(ax, altim_date, geotile_row["$(varname)_altim"] .- mean(geotile_row["$(varname)_altim"][altim_index]), label="Altimetry");
    lines!(ax, dv_altim[geotile=At(single_geotile_test)],  label="Altimetry- RAW")

    axislegend(ax, position=:rt);
    display(f)
end
end