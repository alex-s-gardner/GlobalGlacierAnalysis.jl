using ArchGDAL

folder = "/Users/gardnera/data/GlacierOutlines/RGI2000-v7.0-C-global"

cp(folder, folder * "-fix", force=true)

shp = allfiles(folder; subfolders=true, fn_endswith=".shp")

for s in shp
    d = ArchGDAL.read(s)
    fn_out = replace(s, folder => (folder * "-fix"))
    ArchGDAL.unsafe_gdalvectortranslate([d], ["-nlt" "MULTIPOLYGON" "-dim" "XY"], dest=fn_out)
end