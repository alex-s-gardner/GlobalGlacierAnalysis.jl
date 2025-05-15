
"""
    download_black_marble_2016()

Download NASA's Black Marble 2016 nighttime imagery tiles.

Downloads the 8 tiles (A1-D2) of the Black Marble 2016 nighttime imagery from NASA's
Earth Observatory website. Uses aria2c for efficient downloading with automatic retry
capability in case of connection failures.

The files are saved to the directory "/Users/gardnera/data/BaseMaps/BlackMarble_2016".
"""
using Aria2_jll

finished = false
while finished == false
    try
        for i = 'A':'D', j = 1:2
            fn = "https://eoimages.gsfc.nasa.gov/images/imagerecords/144000/144898/BlackMarble_2016_$(i)$(j)_geo.tif"
            folder = "/Users/gardnera/data/BaseMaps/BlackMarble_2016"

            cmd = `$(Aria2_jll.aria2c()) $fn -c -d $folder`
            run(cmd)
        end
        finished = true
    catch
        println("Failed to download, retrying...")
    end
end