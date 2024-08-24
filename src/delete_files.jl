using Dates
folder = "/mnt/bylot-r3/data/binned/2deg"
#folder = "/mnt/bylot-r3/data/binned_unfiltered/2deg"
f = Altim.allfiles(folder; fn_contains = ["filled"] )

last_modified = Dates.unix2datetime.(mtime.(f))
delete_index = last_modified .< DateTime(2024, 7, 17)

rm.(f[delete_index])