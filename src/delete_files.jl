"""
    delete_files.jl

Delete outdated files from a data directory.

This script removes files containing "filled" in their names that were last modified
before July 17, 2024. It operates on files in the binned 2-degree data directory.
"""

using Dates
import GlobalGlacierAnalysis as GGA

folder = GGA.analysis_paths(; geotile_width).binned

f = GGA.allfiles(folder; fn_contains=["filled"])

last_modified = Dates.unix2datetime.(mtime.(f))
delete_index = last_modified .< DateTime(2024, 7, 17)

rm.(f[delete_index])

