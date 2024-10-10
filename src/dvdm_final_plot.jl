#begin
    using Altim
    using FileIO
    using CairoMakie
    using StatsBase
    using ColorSchemes
    using DataFrames

    geotile_width = 2;
    project_id = :v01
    daterange = (2000, 2024)

    dvdm_synthesis_id = "glacier_filtered"

    final_data_dir = joinpath(Altim.pathlocal.data_dir, "altim_final", "$(geotile_width)_$(project_id)")
    final_filename = "dvdm_$(dvdm_synthesis_id).jld2"
    #final_filename = "dvdm.jld2"
    
    final_figure_dir = joinpath(Altim.pathlocal.data_dir, "altim_final", "$(geotile_width)_$(project_id)", "figures")

    df = load(joinpath(final_data_dir, final_filename), "df")

    df = Altim.dvdm_crop2dates!(df, daterange)

    # recompute trends after crop2dates!(df, daterange)
    df = Altim.geotile_dvdm_add_trend!(df; iterations=1000)

    # rate of acceleration slowing
    regions0 = reduce(vcat, (["rgi$i" for i in vcat(1:12, 16:19)], ["hma"]));
    palette0 = (; color=Makie.resample_cmap(:Dark2_4, length(regions0)))
    
    # all regions together without missions shown
    variable = "dm"
    showmissions = false
    featured_mission = "icesat2"
    (f, regions0, region_offsets, ylims) = Altim.plot_multiregion_dvdm(df;
        variable,
        regions = regions0,
        featured_mission,
        showlines=false,
        showmissions,
        fontsize=15,
        cmap=:Dark2_4,
        title="$(Altim.var2label[variable]): $(Altim.mission2label[featured_mission])"
    )

    begin
        for variable in ["dm", "dv"]
            for featured_mission in unique(df.mission)
            #featured_mission = "grace"

                index = (df.var .==  variable) .& (df.mission .== featured_mission)
                if !any(index)
                    continue
                end
                regions = intersect(unique(df.rgi[index]), regions0)
                index = fill(0, length(regions))
                c = 1;
                for (i,r) in enumerate(regions0)
                    if any(regions .== r)
                        index[c] = i
                        c += 1
                    end
                end
                
                regions = regions0[index]
                palette = (; color=palette0.color[index])

                (f, _, region_offsets, ylims) = Altim.plot_multiregion_dvdm(df;
                    variable,
                    regions,
                    featured_mission,
                    showlines=false,
                    showmissions,
                    fontsize=15,
                    cmap=:Dark2_4,
                    palette,
                    regions_ordered=true,
                    title="$(Altim.var2label[variable]): $(Altim.mission2label[featured_mission])"
                )

                display(f)
                filename = "$(variable)_$(featured_mission)_missions=$(showmissions).png"
                save(joinpath(final_figure_dir, filename), f)
            end
        end
    end

    # firn air content
    begin
        for variable in ["fac", "refreeze"]
        featured_mission = "gemb"
        showmissions = false
        regions_ordered=false
        palette = palette0;
        regions = regions0;

        (f, _, region_offsets, ylims) = Altim.plot_multiregion_dvdm(copy(df);
            variable,
            featured_mission,
            regions,
            showlines=false,
            showmissions,
            fontsize=15,
            regions_ordered=true,
            region_offsets = nothing,
            ylims = nothing,
            cmap=:Dark2_4,
            palette,
            title=Altim.var2label[variable]
        )

        display(f)
        filename0 = "$(variable)_$(featured_mission)_missions=$(showmissions).png"
        save(joinpath(final_figure_dir, filename0), f)
        end
    end


    # smb
    begin
        variable = "smb"
        featured_mission = "gemb"
        showmissions = false
        regions_ordered=true
        regions = regions0;

        (f, _, region_offsets, ylims) = Altim.plot_multiregion_dvdm(copy(df);
            variable,
            featured_mission,
            regions,
            showlines=false,
            showmissions,
            fontsize=15,
            regions_ordered=false,
            region_offsets = nothing,
            ylims = nothing,
            cmap=:Dark2_4,
            title=Altim.var2label[variable],
            delta_offset = -100
        )

        display(f)
        filename = "$(variable)_$(featured_mission)_missions=$(showmissions).png"
        save(joinpath(final_figure_dir, filename), f)
    end

    begin
        variable = "runoff"
        featured_mission = "gemb"
        showmissions = false
        regions_ordered = true
        regions = regions0;

        (f, regions, region_offsets, ylims) = Altim.plot_multiregion_dvdm(copy(df);
            variable,
            featured_mission,
            regions,
            showlines=false,
            showmissions,
            fontsize=15,
            regions_ordered=false,
            region_offsets=nothing,
            ylims=nothing,
            cmap=:Dark2_4,
            title=Altim.var2label[variable],
            delta_offset=-200
        )

        display(f)
        filename = "$(variable)_$(featured_mission)_missions=$(showmissions).png"
        save(joinpath(final_figure_dir, filename), f)
    end

    begin
        variable = "melt"
        featured_mission = "gemb"
        showmissions = false
        regions_ordered = true
 regions = regions0;

        (f, regions, region_offsets, ylims) = Altim.plot_multiregion_dvdm(copy(df);
            variable,
            featured_mission,
            regions,
            showlines=false,
            showmissions,
            fontsize=15,
            regions_ordered=false,
            region_offsets=nothing,
            ylims=nothing,
            cmap=:Dark2_4,
            title=Altim.var2label[variable],
            delta_offset=-200
        )

        display(f)
        filename = "$(variable)_$(featured_mission)_missions=$(showmissions).png"
        save(joinpath(final_figure_dir, filename), f)
    end


    palette0 = (; color=Makie.resample_cmap(:Dark2_4, length(missions)));

    dates = DataFrames.metadata(df, "date")
    decyear = Altim.decimalyear.(dates)

   rgis = unique(df.rgi)


   # missions are listed in the oder that they will 0be plotted
    missions = ["hugonnet", "icesat", "gedi", "icesat2", "gemb","grace", "synthesis"]
    colors = ColorSchemes.tab10;

    for rgi in rgis
    #rgi = "rgi1"
    variable = "dm"
    f = Altim.plot_regional_dvdm(df;
            rgi,
        variable,
        missions,
        colors,
        xticks = 2000:2:2024,
    )
    display(f)
    end


    df0 = Altim.dvdm_bin(df; bin_edges = 2000:.25:2024)

    df0 = Altim.dvdm_delta(df0; anomaly = false)

    df0 = Altim.dvdm_stairs(df0)

    variable = "dm"

    rgis = unique(df.rgi)

    missions = ["synthesis"]; #["hugonnet", "icesat", "gedi", "icesat2", "gemb","grace", "synthesis"]
    colors = ColorSchemes.tab10;
    for rgi in rgis
        f = Altim.plot_regional_dvdm(df0;
                rgi,
            variable,
            missions,
            colors,
            xticks = 2000:2:2024,
        )
        display(f)
    end

    featured_mission = "gemb"
    showmissions = false
    (f, regions0, region_offsets, ylims) = Altim.plot_multiregion_dvdm(df0;
        variable,
        regions = regions0,
        featured_mission,
        showlines=false,
        showmissions,
        fontsize=15,
        cmap=:Dark2_4,
        title="$(Altim.var2label[variable]): $(Altim.mission2label[featured_mission])"
    );
display(f)
end