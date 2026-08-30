using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using Dates
using LsqFit: curve_fit
using Statistics: mean

@testset "Parametric Model Functions" begin
    @testset "model3 - trend + seasonality" begin
        # Model: offset + slope*t + accel*t² + amp*sin(2π*(t+phase))
        t = [0.0, 0.25, 0.5, 0.75]
        params = [5.0, 0.5, 0.0, 2.0, 0.0]  # offset=5, slope=0.5, accel=0, amp=2, phase=0

        result = GGA.model3(t, params)

        # At t=0: 5 + 0 + 0 + 2*sin(0) = 5.0
        @test result[1] ≈ 5.0 rtol=1e-6

        # At t=0.25: 5 + 0.5*0.25 + 0 + 2*sin(2π*0.25) = 5.125 + 2*sin(π/2) = 5.125 + 2 = 7.125
        @test result[2] ≈ 7.125 rtol=1e-6

        # At t=0.5: 5 + 0.5*0.5 + 0 + 2*sin(π) = 5.25 + 0 = 5.25
        @test result[3] ≈ 5.25 rtol=1e-6

        # At t=0.75: 5 + 0.5*0.75 + 0 + 2*sin(3π/2) = 5.375 - 2 = 3.375
        @test result[4] ≈ 3.375 rtol=1e-6
    end

    @testset "model4 - offset + seasonality" begin
        # Model: offset + amp*sin(2π*(t+phase))
        t = [0.0, 0.25, 0.5, 0.75]
        params = [10.0, 3.0, 0.0]  # offset=10, amp=3, phase=0

        result = GGA.model4(t, params)

        # At t=0: 10 + 3*sin(0) = 10.0
        @test result[1] ≈ 10.0 rtol=1e-6

        # At t=0.25: 10 + 3*sin(π/2) = 10 + 3 = 13.0
        @test result[2] ≈ 13.0 rtol=1e-6

        # At t=0.5: 10 + 3*sin(π) = 10.0
        @test result[3] ≈ 10.0 rtol=1e-6

        # At t=0.75: 10 + 3*sin(3π/2) = 10 - 3 = 7.0
        @test result[4] ≈ 7.0 rtol=1e-6
    end

    @testset "offset_trend_seasonal" begin
        # Model: offset + trend*t + amp*sin(2π*(t+phase))
        t = [0.0, 0.5, 1.0]
        params = [1.0, 0.5, 2.0, 0.0]  # offset=1, trend=0.5, amp=2, phase=0

        result = GGA.offset_trend_seasonal(t, params)

        # At t=0: 1 + 0 + 2*sin(0) = 1.0
        @test result[1] ≈ 1.0 rtol=1e-6

        # At t=0.5: 1 + 0.5*0.5 + 2*sin(π) = 1.25 + 0 = 1.25
        @test result[2] ≈ 1.25 rtol=1e-6

        # At t=1.0: 1 + 0.5*1.0 + 2*sin(2π) = 1.5 + 0 = 1.5
        @test result[3] ≈ 1.5 rtol=1e-6
    end

    @testset "offset_trend_seasonal2 - cos/sin form" begin
        # Model: offset + trend*t + c*cos(2πt) + s*sin(2πt)
        t = [0.0, 0.25, 0.5, 0.75]
        params = [5.0, 1.0, 2.0, 1.0]  # offset=5, trend=1, cos_coef=2, sin_coef=1

        result = GGA.offset_trend_seasonal2(t, params)

        # At t=0: 5 + 0 + 2*cos(0) + 1*sin(0) = 5 + 2 + 0 = 7.0
        @test result[1] ≈ 7.0 rtol=1e-6

        # At t=0.25: 5 + 0.25 + 2*cos(π/2) + 1*sin(π/2) = 5.25 + 0 + 1 = 6.25
        @test result[2] ≈ 6.25 rtol=1e-6

        # At t=0.5: 5 + 0.5 + 2*cos(π) + 1*sin(π) = 5.5 - 2 + 0 = 3.5
        @test result[3] ≈ 3.5 rtol=1e-6
    end

    @testset "model1_trend - linear elevation model" begin
        # Model: intercept + slope*h
        h = [0.0, 1000.0, 2000.0]
        params = [10.0, -0.005]  # intercept=10, slope=-0.005

        result = GGA.model1_trend(h, params)

        # At h=0: 10 + 0 = 10.0
        @test result[1] ≈ 10.0 rtol=1e-6

        # At h=1000: 10 - 0.005*1000 = 10 - 5 = 5.0
        @test result[2] ≈ 5.0 rtol=1e-6

        # At h=2000: 10 - 0.005*2000 = 10 - 10 = 0.0
        @test result[3] ≈ 0.0 rtol=1e-6
    end

    @testset "model2 - quadratic elevation model" begin
        # Model: p[1] + p[2]*h + p[3]*h²
        h = [0.0, 1000.0, 2000.0]
        params = [100.0, -0.05, 0.00001]  # intercept, linear, quadratic

        result = GGA.model2(h, params)

        # At h=0: 100
        @test result[1] ≈ 100.0 rtol=1e-6

        # At h=1000: 100 - 0.05*1000 + 0.00001*1000^2 = 100 - 50 + 10 = 60.0
        @test result[2] ≈ 60.0 rtol=1e-5

        # At h=2000: 100 - 0.05*2000 + 0.00001*2000^2 = 100 - 100 + 40 = 40.0
        @test result[3] ≈ 40.0 rtol=1e-5
    end

    @testset "offset_trend_acceleration_seasonal2" begin
        # Model: offset + linear*t + quad*t² + cos_coef*cos(2πt) + sin_coef*sin(2πt)
        t = [0.0, 0.5, 1.0]
        params = [10.0, 2.0, 0.5, 1.0, 0.5]  # offset, linear, quad, cos, sin

        result = GGA.offset_trend_acceleration_seasonal2(t, params)

        # At t=0: 10 + 0 + 0 + 1*1 + 0.5*0 = 11.0
        @test result[1] ≈ 11.0 rtol=1e-6

        # At t=0.5: 10 + 2*0.5 + 0.5*0.25 + 1*cos(π) + 0.5*sin(π)
        #          = 10 + 1 + 0.125 - 1 + 0 = 10.125
        @test result[2] ≈ 10.125 rtol=1e-6

        # At t=1.0: 10 + 2*1 + 0.5*1 + 1*cos(2π) + 0.5*sin(2π)
        #          = 10 + 2 + 0.5 + 1 + 0 = 13.5
        @test result[3] ≈ 13.5 rtol=1e-6
    end

    @testset "seasonal_peak_fraction" begin
        # A pure cosine peaks at the start of the year, a pure sine a quarter year later.
        @test GGA.seasonal_peak_fraction(1.0, 0.0) ≈ 0.0 atol=1e-12
        @test GGA.seasonal_peak_fraction(0.0, 1.0) ≈ 0.25 atol=1e-12
        @test GGA.seasonal_peak_fraction(-1.0, 0.0) ≈ 0.5 atol=1e-12
        @test GGA.seasonal_peak_fraction(0.0, -1.0) ≈ 0.75 atol=1e-12

        # Always in [0, 1)
        @test all(0 <= GGA.seasonal_peak_fraction(cos(θ), sin(θ)) < 1 for θ in range(-2π, 2π, length=97))

        # Amplitude does not affect the peak timing
        @test GGA.seasonal_peak_fraction(3.0, 3.0) ≈ GGA.seasonal_peak_fraction(0.5, 0.5) atol=1e-12

        # Recover a peak planted at a known day of year by fitting offset_trend_seasonal2.
        dates = DateTime(2010, 1, 15):Month(1):DateTime(2015, 12, 15)
        t = GGA.decimalyear.(collect(dates))
        for peak_day in (15, 100, 200, 300)
            y = 3.0 .* cos.(2π .* (t .- peak_day / 365.25))
            fit = curve_fit(GGA.offset_trend_seasonal2, t .- ceil(mean(t)), y, zeros(4))
            recovered = 365.25 * GGA.seasonal_peak_fraction(fit.param[3], fit.param[4])
            @test recovered ≈ peak_day rtol=1e-4
        end
    end

    @testset "Model vectorization" begin
        # Test that models handle vector inputs correctly
        t_vec = range(0, 1, length=100)
        params_simple = [0.0, 1.0, 0.0, 0.0, 0.0]

        result = GGA.model3(t_vec, params_simple)
        @test length(result) == 100
        @test result[1] ≈ 0.0 rtol=1e-6
        @test result[end] ≈ 1.0 rtol=1e-6
    end
end
