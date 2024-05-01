using Downloads, NetCDF, NCDatasets
path2file = "https://its-live-data.s3.amazonaws.com/Test/N52W175.nc";
Downloads.download(path2file, "N52W175.nc");

# using NetCDF
@time nc = copy(ncread("N52W175.nc", "z"));

# using NCDatasets
@time nc2 = copy(Dataset("N52W175.nc", "r")["z"]);