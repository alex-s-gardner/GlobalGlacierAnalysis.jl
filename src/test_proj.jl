


obj = Proj.proj_create("EPSG:4326")

i1 = C_NULL
i2 = C_NULL

a = Proj.proj_identify(obj, i1, i2)

