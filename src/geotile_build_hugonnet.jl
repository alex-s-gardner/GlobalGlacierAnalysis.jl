# build geotiles from Hugonnet stacks
using Altim
using Extents

force_remake = true; 
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
    
# build geotiles
begin
    printstyled("building Hugonnet geotiles\n"; color=:blue, bold=true)
    for geotile in eachrow(geotiles)
        geotile_build_hugonnet(geotile, paths.hugonnet.geotile, hstacks; force_remake=force_remake)
    end
end
