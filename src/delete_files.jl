"""
Delete old data files from a specified folder based on modification date.

This script removes files that were last modified before a cutoff date.
It specifically targets files containing "filled" in their names.

The script:
1. Finds all files in the target folder containing "filled" 
2. Gets last modified timestamps for each file
3. Creates boolean mask for files older than cutoff date
4. Deletes files matching the date criteria

Parameters:
- folder: Path to folder containing files to check
- cutoff_date: Files modified before this date will be deleted (2024-07-17)

Dependencies:
Dates, Altim
"""

using Dates
folder = "/mnt/bylot-r3/data/binned/2deg"
#folder = "/mnt/bylot-r3/data/binned_unfiltered/2deg"
f = Altim.allfiles(folder; fn_contains = ["filled"] )

last_modified = Dates.unix2datetime.(mtime.(f))
delete_index = last_modified .< DateTime(2024, 7, 17)

rm.(f[delete_index])