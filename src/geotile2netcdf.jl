


rgi = [2, 3, 4]
reference_run = "/mnt/bylot-r3/data/binned_unfiltered/2deg/glacier_dh_best_cc_meanmadnorm3_v01_filled_ac_p1_synthesized.jld2"
reference_run_nc = replace(reference_run, ".jld2" => ".nc")


#for binned_synthesized_file in path2runs
binned_synthesized_file = reference_run

binned_synthesized_dv_file = replace(binned_synthesized_file, ".jld2" => "_gembfit_dv.jld2")
geotiles0 = FileIO.load(binned_synthesized_dv_file, "geotiles")


#subset to rgi
df = geotiles0[in.(geotiles0.rgiid, Ref(rgi)), :]
valid_model_data = .!isnan.(df[1, :runoff])
valid_altim_data = .!isnan.(df[1, :dm_altim])

# package as netcdf
NCDataset(reference_run_nc, "c") do ds

    # First all dimensions need to be defined
    defDim(ds, "geotile", nrow(df))
    defDim(ds, "time_model", sum(valid_model_data))
    defDim(ds, "time_altim", sum(valid_altim_data))

    # now add the variables
    v = defVar(ds, "geotile", df.id, ("geotile",))
    v.attrib["long_name"] = "geotile id"
    v.attrib["notes"] = "defines the extents of the geotile"

    v = defVar(ds, "time_model", colmetadata(df, :runoff, "date")[valid_model_data], ("time_model",))
    v.attrib["long_name"] = "time"
    v.attrib["notes"] = "time of model data"

    v = defVar(ds, "time_altim", colmetadata(df, :dv_altim, "date")[valid_altim_data], ("time_altim",))
    v.attrib["long_name"] = "time"
    v.attrib["notes"] = "time of altimetry data"

    # 1d variables
    v = defVar(ds, "area", sum.(df.area_km2), ("geotile",))
    v.attrib["long_name"] = "glacier area"
    v.attrib["units"] = "km2"

    # 2d variables
    vars = ["runoff", "fac", "smb", "rain", "acc", "melt", "ec", "refreeze", "dv", "dm", "dm_altim", "dv_altim"]
    for var in vars 
        if occursin("altim", var)
            data = reduce(hcat, df[:, var])[valid_altim_data, :]'
            tvar = "time_altim"
            long_name = "$var: from altimetry synthesis"
        else
            data = reduce(hcat, df[:, var])[valid_model_data, :]'
            tvar = "time_model"
            long_name = "$var: from calibrated model output"
        end

        v = defVar(ds, var, data, ("geotile", tvar))
        v.attrib["units"] = "m.w.e."
        v.attrib["long_name"] = long_name
    end
end