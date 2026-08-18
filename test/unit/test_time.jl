using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
using Dates

@testset "Time Conversions" begin
    @testset "decimalyear basic conversion" begin
        # Test year boundary
        dt_jan1 = DateTime(2019, 1, 1, 0, 0, 0)
        @test GGA.decimalyear(dt_jan1) == 2019.0

        # Test mid-year (non-leap year)
        # July 1 is day 182 of 365
        dt_july1 = DateTime(2019, 7, 1)
        expected = 2019.0 + (182 / 365)
        @test GGA.decimalyear(dt_july1) ≈ expected rtol=1e-6

        # Test end of year
        dt_dec31 = DateTime(2019, 12, 31, 23, 59, 59)
        @test GGA.decimalyear(dt_dec31) ≈ 2020.0 rtol=1e-4
    end

    @testset "decimalyear leap year" begin
        # Leap year has 366 days
        dt_july1_leap = DateTime(2020, 7, 1)
        # July 1 is day 183 of 366 in a leap year
        expected = 2020.0 + (183 / 366)
        @test GGA.decimalyear(dt_july1_leap) ≈ expected rtol=1e-6

        # End of leap year
        dt_dec31_leap = DateTime(2020, 12, 31)
        expected_end = 2020.0 + (366 / 366)
        @test GGA.decimalyear(dt_dec31_leap) ≈ expected_end rtol=1e-4
    end

    @testset "decimalyear2datetime basic conversion" begin
        # Test year boundary
        dt = GGA.decimalyear2datetime(2019.0)
        @test Dates.year(dt) == 2019
        @test Dates.month(dt) == 1
        @test Dates.day(dt) == 1

        # Test mid-year
        dt_mid = GGA.decimalyear2datetime(2019.5)
        @test Dates.year(dt_mid) == 2019
        # Mid-year should be around July 2 (day 183 of 365)
        @test Dates.month(dt_mid) >= 6 && Dates.month(dt_mid) <= 7
    end

    @testset "decimalyear2datetime leap year" begin
        # Test leap year conversion
        dt_leap = GGA.decimalyear2datetime(2020.5)
        @test Dates.year(dt_leap) == 2020
        @test Dates.isleapyear(dt_leap) == true
        # Mid-year in leap year should also be around July
        @test Dates.month(dt_leap) >= 6 && Dates.month(dt_leap) <= 7
    end

    @testset "Round-trip conversion" begin
        # Test round-trip for various dates
        test_dates = [
            DateTime(2018, 1, 1),
            DateTime(2018, 6, 15),
            DateTime(2019, 12, 31),
            DateTime(2020, 2, 29),  # Leap day
            DateTime(2021, 7, 4),
            DateTime(2022, 10, 24, 15, 30, 45)
        ]

        for dt in test_dates
            decyear = GGA.decimalyear(dt)
            dt_back = GGA.decimalyear2datetime(decyear)

            # Allow small tolerance due to floating point and millisecond rounding
            diff_seconds = abs(Dates.value(dt - dt_back)) / 1000
            @test diff_seconds < 1.0  # Within 1 second
        end
    end

    @testset "Round-trip many random dates" begin
        # Test 100 random dates over 10 years
        Random.seed!(42)
        start_date = DateTime(2015, 1, 1)
        end_date = DateTime(2025, 1, 1)
        date_range_ms = Dates.value(end_date - start_date)

        for _ in 1:100
            # Generate random DateTime
            random_ms = rand(0:date_range_ms)
            dt = start_date + Millisecond(random_ms)

            # Round-trip
            decyear = GGA.decimalyear(dt)
            dt_back = GGA.decimalyear2datetime(decyear)

            # Check within 1 second tolerance
            diff_seconds = abs(Dates.value(dt - dt_back)) / 1000
            @test diff_seconds < 1.0
        end
    end

    @testset "Edge cases" begin
        # Test century boundary
        dt_2000 = DateTime(2000, 1, 1)
        @test GGA.decimalyear(dt_2000) == 2000.0

        # Test far future
        dt_future = DateTime(2100, 1, 1)
        @test GGA.decimalyear(dt_future) == 2100.0

        # Test subsecond precision
        dt_precise = DateTime(2020, 1, 1, 12, 30, 45, 123)
        decyear = GGA.decimalyear(dt_precise)
        @test decyear > 2020.0
        @test decyear < 2020.1
    end
end
