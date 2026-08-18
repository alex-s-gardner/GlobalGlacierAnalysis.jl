"""
Synthetic river network generation for testing routing functions
"""

using DataFrames

"""
    create_linear_network(n_nodes=5)

Create a simple linear river network: 1 → 2 → 3 → ... → n → 0 (ocean)
"""
function create_linear_network(n_nodes=5)
    ids = collect(1:n_nodes)
    next_ids = [i < n_nodes ? i+1 : 0 for i in ids]
    lengths_km = fill(10.0, n_nodes)
    lons = range(-120, -110, length=n_nodes)
    lats = fill(45.0, n_nodes)

    return DataFrame(
        COMID=ids,
        NextDownID=next_ids,
        lengthkm=lengths_km,
        lon=lons,
        lat=lats
    )
end

"""
    create_branching_network()

Create a branching network:
    1 → 3 ← 2
        ↓
        4 → 5 → 0
"""
function create_branching_network()
    ids = [1, 2, 3, 4, 5]
    next_ids = [3, 3, 4, 5, 0]
    lengths_km = [10.0, 10.0, 15.0, 20.0, 25.0]
    lons = [-120.0, -119.0, -118.0, -116.0, -114.0]
    lats = [45.5, 44.5, 45.0, 45.0, 45.0]

    return DataFrame(
        COMID=ids,
        NextDownID=next_ids,
        lengthkm=lengths_km,
        lon=lons,
        lat=lats
    )
end

"""
    create_dendritic_network()

Create a more complex dendritic (tree-like) network:
    1 → 4 ← 2
        ↓
    3 → 5 → 6 → 0
"""
function create_dendritic_network()
    ids = [1, 2, 3, 4, 5, 6]
    next_ids = [4, 4, 5, 5, 6, 0]
    lengths_km = [8.0, 12.0, 10.0, 15.0, 20.0, 30.0]
    lons = [-120.0, -119.5, -118.5, -119.0, -117.0, -115.0]
    lats = [46.0, 45.5, 44.5, 45.5, 45.0, 45.0]

    return DataFrame(
        COMID=ids,
        NextDownID=next_ids,
        lengthkm=lengths_km,
        lon=lons,
        lat=lats
    )
end
