"""
Download NASA's Black Marble nighttime lights imagery from 2016.

This script downloads the Black Marble nighttime lights imagery tiles from NASA's servers.
The imagery is split into 8 tiles (A1-D2) covering different regions of the Earth.

The script:
1. Uses aria2c for reliable downloads with resume capability
2. Retries failed downloads automatically
3. Downloads all 8 tiles (A1-D2) to the specified local folder

Dependencies:
- Aria2_jll: For the aria2c download utility

Output:
- Downloads .tif files to the specified folder path
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