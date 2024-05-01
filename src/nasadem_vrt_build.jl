# build NASADEM vrt
using ArchGDAL, Altim

path2folder = "/mnt/devon-r2/shared_data/NASADEM/mosaic/";
suffix = ["wgs84_hgt.tif", "egm08_hgt.tif", "num.tif", "swb.tif"];

# build single .vrt for all .tif files with suffix
for suff in suffix
    suff = first(suffix)
    out_vrt = joinpath(path2folder, replace(suff, ".tif" => ".vrt"))
    in_tifs =  searchdir(path2folder, suff);
    in_tifs = joinpath.(path2folder, in_tifs)
    in_tifs = ArchGDAL.read.(in_tifs)

    ArchGDAL.gdalbuildvrt(in_tifs; dest=out_vrt) do vrt
        ArchGDAL.write(vrt, out_vrt)
    end
end

