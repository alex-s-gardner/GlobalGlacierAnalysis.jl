# build geotiles from Hugonnet stacks
using Altim
using Extents

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

#adding x/y offset to Hugonnet geotiles 

# offsets added to all hugonnet data on 2023/08/21
mission = :hugonnet
xoffset = +50;
yoffset = -50;
Altim.geotile_addxy_offsets(geotiles, paths[mission].geotile, xoffset, yoffset)