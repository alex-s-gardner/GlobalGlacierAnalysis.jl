using ModelParameters
using Distributions

@kwdef struct MyStruct00
    a = Param([0,1,3]; )

end


using ModelParameters

Base.@kwdef struct Project{A}
    update_geotile_missions::A = Param(["icesat", "icesat2"],)
end

model = Model(Project(), )


Base.@kwdef struct Project{A,B}
    project_id::A = Param(0.8, bounds=(0.2, 0.9))
    geotile_width::B = Param(0.5, bounds=(0.7, 0.4))
end

Base.@kwdef struct Binning{A, B, C}  
    force_remake_binning::A =  Param(false,);
    update_geotile::B = Param(false,); # this will load in prevous results to update select geotiles or missions
    update_geotile_missions::C = Param(["icesat", "icesat2"],)
end

model = Model((Project(), Binning()))


# run parameters
all_permutations_for_glacier_only = true
surface_masks = [:glacier, :glacier_rgi7, :land, :glacier_b1km, :glacier_b10km]
hugonnet_filtered_flags = [true, false]
dem_ids = [:best, :cop30_v2]
binning_methods = ["meanmadnorm3", "meanmadnorm5", "median", "meanmadnorm10"]
curvature_corrects = [true, false]
max_canopy_height = 1
filt = (;dh_max=200)






Base.@kwdef struct Submodel1{A,B}
    ::A = Param(0.8, bounds=(0.2, 0.9))
    β::B = Param(0.5, bounds=(0.7, 0.4))
end

Base.@kwdef struct Submodel2{Γ}
    γ::Γ = Param(1e-3, bounds=(1e-4, 1e-2))
end

Base.@kwdef struct SubModel3{Λ,X}
    λ::Λ = Param(0.8, bounds=(0.2, 0.9))
    x::X = Submodel2()
end

model = Model((Submodel1(), SubModel3()))



Base.@kwdef struct GeotileBinning{A,B}
    project_id::A = Param(0.8, bounds=(0.2, 0.9))
    geotile_width::B = Param(0.5, bounds=(0.7, 0.4))
end

model = ModelParameters.Model((GeotileBinning,GeotileBinning))


 project_id = :v01
    geotile_width = 2

    warnings = false
    showplots = false

    # run parameters
    force_remake_binning = false;
    update_geotile = false; # this will load in prevous results to update select geotiles or missions
    update_geotile_missions = ["icesat2"]

    # run parameters
    all_permutations_for_glacier_only = true
    surface_masks = [:glacier, :glacier_rgi7, :land, :glacier_b1km, :glacier_b10km]
    hugonnet_filtered_flags = [true, false]
    dem_ids = [:best, :cop30_v2]
    binning_methods = ["meanmadnorm3", "meanmadnorm5", "median", "meanmadnorm10"]
    curvature_corrects = [true, false]
    max_canopy_height = 1
    filt = (;dh_max=200)

    # filling only parameters
    filling_paramater_sets = [1, 2]
    amplitude_corrects = [true, false]
    force_remake_fill  = false
    plot_dh_as_function_of_time_and_elevation = true;
    mission_reference_for_amplitude_normalization = "icesat2"



 project_id,
        geotile_width,
        warnings,
        showplots,

        # run parameters
        force_remake = force_remake_binning,
        update_geotile, # this will load in prevous results to update select geotiles or missions
        update_geotile_missions,

        # run parameters
        all_permutations_for_glacier_only,
        surface_masks,
        hugonnet_filtered_flags,
        dem_ids,
        binning_methods,
        curvature_corrects,

        #### DON NOT CHANGE THESE PARAMETERS
        max_canopy_height, # do not change

        # filter parameters
        filt,
        

end