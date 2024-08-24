# consider running `geotile_remove_files.jl` to clear files from folder... if there is an 
# error this will allow you to restart processing where you left off, unlike 
# force_remake = true


include("/home/gardnera/Documents/GitHub/Altim.jl/src/geotile_build_archive.jl")
include("/home/gardnera/Documents/GitHub/Altim.jl/src/geotile_build_hugonnet.jl")
include("/home/gardnera/Documents/GitHub/Altim.jl/src/geotile_dem_extract.jl")
include("/home/gardnera/Documents/GitHub/Altim.jl/src/geotile_mask_extract.jl")
include("/home/gardnera/Documents/GitHub/Altim.jl/src/geotile_canopyh_extract.jl")
include("/home/gardnera/Documents/GitHub/Altim.jl/src/geotile_hyps.jl")


include("/home/gardnera/Documents/GitHub/Altim.jl/src/gemb_classes_binning.jl")
include("/home/gardnera/Documents/GitHub/Altim.jl/src/geotile_binning.jl")
