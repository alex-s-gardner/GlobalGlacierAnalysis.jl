"""
Synthetic geotile data generation for testing
"""

using DimensionalData
import DimensionalData as DD
using Dates

"""
    create_synthetic_geotile_data(; n_dates=12, n_heights=10, trend=-0.5, seasonal_amp=0.2)

Generate synthetic geotile data with known trend and seasonality.

# Returns
- DimArray with dimensions (date, height) and known temporal/elevation patterns
"""
function create_synthetic_geotile_data(; n_dates=12, n_heights=10, trend=-0.5, seasonal_amp=0.2)
    dates = [DateTime(2018,1,1) + Month(i) for i in 0:n_dates-1]
    heights = range(1000, 3000, length=n_heights)

    data = zeros(n_dates, n_heights)

    for (i, date) in enumerate(dates)
        t = (Dates.value(date - DateTime(2018,1,1))) / (365.25 * 24 * 60 * 60 * 1000)  # Years

        for (j, h) in enumerate(heights)
            # Trend + seasonality + elevation dependency + noise
            signal = trend * t
            signal += seasonal_amp * sin(2π * t)
            signal += (h - 2000) * 0.0001  # Slight elevation gradient
            signal += randn() * 0.05  # Noise

            data[i, j] = signal
        end
    end

    return DimArray(
        data,
        (
            DD.Dim{:date}(dates),
            DD.Dim{:height}(heights)
        )
    )
end

"""
    create_synthetic_dem(; size=(100, 100), min_elev=0, max_elev=3000)

Generate synthetic DEM with realistic terrain.
"""
function create_synthetic_dem(; size=(100, 100), min_elev=0, max_elev=3000)
    # Create simple synthetic terrain (cone shape)
    center_x, center_y = size[1] ÷ 2, size[2] ÷ 2

    dem = zeros(size...)

    for i in 1:size[1]
        for j in 1:size[2]
            dist = sqrt((i - center_x)^2 + (j - center_y)^2)
            max_dist = sqrt(center_x^2 + center_y^2)

            # Elevation decreases with distance from center
            dem[i, j] = max_elev * (1 - dist / max_dist) + min_elev
        end
    end

    return dem
end
