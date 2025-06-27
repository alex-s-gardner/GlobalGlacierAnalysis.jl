import GeoFormatTypes as GFT
using ArchGDAL
using LibGEOS
using Proj



proj_string = GFT.ProjString("+proj=utm +zone=1 +datum=WGS84 +units=m +no_defs +type=crs")
crs = Proj.identify(Proj.convert(Proj.CRS, proj_string), auth_name="EPSG")
epsg = convert(GFT.EPSG, first(crs).crs)


convert(GFT.EPSG, foo)





GFT.EPSG(proj_string)

GFT.convert(GFT.EPSG, Proj.convert(Proj.CRS, proj_string))

Proj.to

GFT.convert(GFT.EPSG, proj_string)
