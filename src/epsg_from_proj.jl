proj_string = GFT.ProjString("+proj=utm +zone=1 +datum=WGS84 +units=m +no_defs +type=crs")
crs = Proj.identify(Proj.convert(Proj.CRS, proj_string), auth_name="EPSG")
epsg = convert(GFT.EPSG, first(crs).crs)