# repair_gedi_placeholders.jl
#
# Give the GEDI archive's placeholder rows their granule ids back.
#
# Usage:
#   julia --project src/repair_gedi_placeholders.jl --dry-run   # report, change nothing
#   julia --project src/repair_gedi_placeholders.jl             # repair
#
# What is wrong
#
# `geotile_build` records a granule that returned no points inside a geotile by pushing one all-NaN
# placeholder row carrying that granule's id, so the next pass knows not to ask for it again. The code
# that built this archive pushed the placeholder without setting the id, so every one of them reads
# `id == "0"`. 43,846 granules across 267 geotiles are therefore invisible to the incremental rule and
# would be requested again on every future update, forever.
#
# Why the ids are recoverable
#
# `granules.remote` on disk is still the list the archive was built against, so for each geotile
#
#     missing = (granules in the remote list) - (granule ids present in the geotile)
#
# is exactly the set the placeholders stand for. That the two counts agree is checkable per geotile,
# and they do: in 45 of 46 geotiles sampled before writing this, `count(id == "0")` equalled
# `length(missing)` exactly. A geotile where they disagree is left alone and reported -- the update
# will simply re-request those granules, which is correct if slower.
#
# This must run BEFORE the CMR search is re-run. A new search overwrites `granules.remote` with
# granules acquired since, at which point `missing` no longer means "placeholders" and the
# reconstruction is gone. The script snapshots the list and refuses to proceed if a snapshot already
# disagrees with what is on disk.

using Arrow
using DataFrames
using Printf
using ProgressMeter
import GlobalGlacierAnalysis as GGA

const PLACEHOLDER_ID = "0"
const GRANULE_PREFIX = "GEDI02_A_"

dry_run = "--dry-run" in ARGS

project_id = :v01
geotile_width = 2

paths = GGA.project_paths(; project_id)
geotile_dir = paths.gedi.geotile
remote_file = paths.gedi.granules_remote
snapshot_file = remote_file * ".pre-update"

printstyled("GEDI placeholder repair$(dry_run ? " [DRY RUN]" : "")\n"; color=:blue, bold=true)
println("  geotiles: ", geotile_dir)
println("  granule list: ", remote_file)

isfile(remote_file) || error("no remote granule list at $remote_file")

# Snapshot the list the reconstruction depends on, so the repair stays auditable after the search runs.
if isfile(snapshot_file)
    if read(snapshot_file) != read(remote_file)
        error("$snapshot_file already exists and differs from $remote_file -- the granule list has " *
              "been re-searched since the snapshot was taken, so the placeholder mapping can no " *
              "longer be reconstructed. Restore the snapshot over the current list, or skip the " *
              "repair and let the update re-request those granules.")
    end
    println("  snapshot: already present and identical")
elseif dry_run
    println("  snapshot: would write ", snapshot_file)
else
    cp(remote_file, snapshot_file)
    println("  snapshot: wrote ", snapshot_file)
end

remote = Dict(row.id => Set(String(g.id) for g in row.granules)
              for row in eachrow(DataFrame(Arrow.Table(remote_file))))

files = sort(filter(f -> endswith(f, ".arrow"), readdir(geotile_dir)))
println("  geotile files: ", length(files))

"""
    repair_geotile(path, geotile, remote_ids; dry_run=false) -> (repaired_rows, skip_reason)

Re-key one geotile's placeholder rows. Returns the number of rows repaired and, when nothing was
changed, why: `nothing` means the file was already clean.
"""
function repair_geotile(path, geotile, remote_ids; dry_run=false)
    df = try
        DataFrame(Arrow.Table(path))
    catch e
        return (0, (reason="unreadable [$(sprint(showerror, e))]", placeholders=-1, missing_n=-1))
    end

    placeholders = findall(==(PLACEHOLDER_ID), df.id)
    isempty(placeholders) && return (0, nothing)

    present = Set(id for id in unique(df.id) if startswith(id, GRANULE_PREFIX))
    missing_ids = sort(collect(setdiff(remote_ids, present)))

    # The check that makes this safe. Only rewrite when the placeholder count is exactly the number of
    # granules the archive cannot account for -- then the correspondence is forced, not guessed.
    if length(placeholders) != length(missing_ids)
        return (0, (reason="count mismatch", placeholders=length(placeholders),
            missing_n=length(missing_ids)))
    end

    dry_run && return (length(placeholders), nothing)

    ids = copy(df.id)
    ids[placeholders] = missing_ids
    df.id = ids

    GGA.atomic_write(path; suffix=".arrow") do tmp
        Arrow.write(tmp, df)
    end
    return (length(placeholders), nothing)
end

repaired_tiles = 0
repaired_rows = 0
clean_tiles = 0
skipped = NamedTuple[]

progress = Progress(length(files); dt=1, desc="Repairing GEDI placeholders...")

for file in files
    next!(progress)
    geotile = replace(file, ".arrow" => "")

    if !haskey(remote, geotile)
        push!(skipped, (id=geotile, reason="not in the remote granule list",
            placeholders=-1, missing_n=-1))
        continue
    end

    rows, skip = repair_geotile(joinpath(geotile_dir, file), geotile, remote[geotile]; dry_run)

    if !isnothing(skip)
        push!(skipped, (id=geotile, skip...))
    elseif rows == 0
        global clean_tiles += 1
    else
        global repaired_tiles += 1
        global repaired_rows += rows
    end
end
finish!(progress)

@printf("\n%s %d placeholder rows across %d geotiles\n",
    dry_run ? "would repair" : "repaired", repaired_rows, repaired_tiles)
@printf("%d geotiles had no placeholders to repair\n", clean_tiles)

if !isempty(skipped)
    printstyled("\n$(length(skipped)) geotile(s) left untouched:\n"; color=:light_yellow)
    for s in skipped
        if s.placeholders >= 0
            @printf("  %-28s %s (%d placeholders, %d unaccounted granules)\n",
                s.id, s.reason, s.placeholders, s.missing_n)
        else
            @printf("  %-28s %s\n", s.id, s.reason)
        end
    end
    println("\nThese are not errors: the update will re-request the granules involved, which costs " *
            "time but produces the same archive.")
end

if dry_run
    printstyled("\nDry run -- nothing was written. Re-run without --dry-run to apply.\n"; color=:blue)
else
    printstyled("\nRepair complete. The CMR search can now be re-run safely.\n"; color=:green)
end
