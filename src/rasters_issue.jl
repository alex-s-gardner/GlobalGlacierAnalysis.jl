out_folder = "/mnt/bylot-r3/data/junk"

fname = "https://its-live-data.s3.amazonaws.com/velocity_mosaic/v2/static/cog/ITS_LIVE_velocity_120m_RGI01A_0000_v02_v.tif"
rasterfile = joinpath(out_folder, splitpath(fname)[end])

if !isfile(rasterfil)
    Downloads.download(fname, rasterfil)
end


