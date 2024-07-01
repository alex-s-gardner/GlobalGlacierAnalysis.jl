# build geotiles from Hugonnet stacks
using Altim
using Extents

force_remake = false; 
geotile_width = 2; #geotile width [degrees]
paths = project_paths();
geotiles = project_geotiles(; geotile_width=geotile_width);

#= --------------------------------------------------------------------------
@warn "--- only processing a subset of geotiles ----"
region = :WNA;
extent, epsg = region_extent(region);
geotiles = geotile_subset!(geotiles, extent);
=# 

# custome subset
# geotiles = geotile_subset!(geotiles, Extent(X=(-180, 180), Y = (85,90)));

# make directly if it doesn't exist
if !isdir(paths.hugonnet.geotile)
    mkpath(paths.hugonnet.geotile)
end

# get hstack catalogue
hstacks = hstack_catalogue(paths.hugonnet.raw_data; update_catalogue=force_remake)

if contains(hstacks.path[1], "prefilt")
    old_format = false
else
    old_format = true
end

# build geotiles
begin
    printstyled("building Hugonnet geotiles\n"; color=:blue, bold=true)
    for geotile in eachrow(geotiles)
        geotile_build_hugonnet(geotile, paths.hugonnet.geotile, hstacks; force_remake, old_format)
    end
end
