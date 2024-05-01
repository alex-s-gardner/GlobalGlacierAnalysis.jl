using Altim

dataset = :arcticdem_v4_10m
dataset = :rema_v2_10m

baffin_path = ["/","mnt", "baffin-r1"];

oldpath = splitpath(p[dataset]);
newpath = joinpath(vcat(baffin_path,oldpath[4:end-1]))
oldpath = joinpath(oldpath[1:end-1])

mkpath(newpath)
@time Base.Filesystem.cptree(oldpath, newpath, force=true)

