using Altim
p = Altim.project_paths()
folder = p.icesat2.geotile


extent = Extent(Y=(64.,65.), X=(-15.,-13));
suffix = "arrow"
df = Altim.GeoTiles.listtiles(folder; suffix, extent)

d = Altim.GeoTiles.readall(folder; suffix="arrow", extent)
