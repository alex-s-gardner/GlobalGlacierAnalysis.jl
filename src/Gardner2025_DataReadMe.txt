############################################################################################
#                           Gardner2025_geotiles_rates_km3yr.gpkg                          #
#                   m.i.e to m.w.e conversion = 910 kg/m3 / 1000 kg/m3 = 0.91              #
############################################################################################
Variables:
- id: geotile id
- area_km2: glacier area (RGI6) 
- group: geotile grouping id
- pscale: optimal precipitation scaling [unitless]
- Δheight: optimal height offset [m]
- rgiid: RGI region id (RGI6) 
- discharge_trend: solid ice discharge through flux gate located near grounding line [m.i.e./yr]
- dh_synthesis_trend: trend in height anomaly from synthesis of satellite data [m/yr]
- dh_synthesis_amplitude: seasonal amplitude in height anomaly from synthesis of satellite data [m]
- dh_synthesis_phase: day of seasonal maximum in height anomaly from synthesis of satellite data [day of year]
- runoff_trend: trend in glacier runoff [m.i.e./yr]
- runoff_amplitude: seasonal amplitude in glacier runoff [m.i.e.]
- runoff_phase: day of seasonal maximum in glacier runoff [day of year]
- fac_trend: trend in firn air content [m/yr]
- fac_amplitude: seasonal amplitude in firn air content [m]
- fac_phase: day of seasonal maximum in firn air content [day of year]
- smb_trend: trend in surface mass balance [m.i.e./yr]
- smb_amplitude: seasonal amplitude in surface mass balance [m.i.e.]
- smb_phase: day of seasonal maximum in surface mass balance [day of year]
- rain_trend: trend in rainfall [m.i.e./yr]
- rain_amplitude: seasonal amplitude in rainfall [m.i.e.]
- rain_phase: day of seasonal maximum in rainfall [day of year]
- acc_trend: trend in accumulation [m.i.e./yr]
- acc_amplitude: seasonal amplitude in accumulation [m.i.e.]
- acc_phase: day of seasonal maximum in accumulation [day of year]
- melt_trend: trend in surface melt [m.i.e./yr]
- melt_amplitude: seasonal amplitude in surface melt [m.i.e.]
- melt_phase: day of seasonal maximum in surface melt [day of year]
- ec_trend: trend in evaporation/condensation [m.i.e./yr]
- ec_amplitude: seasonal amplitude in evaporation/condensation [m.i.e.]
- ec_phase: day of seasonal maximum in evaporation/condensation [day of year]
- refreeze_trend: trend in refreezing [m.i.e./yr]
- refreeze_amplitude: seasonal amplitude in refreezing [m.i.e.]
- refreeze_phase: day of seasonal maximum in refreezing [day of year]
- dh_trend: trend in height change [m/yr]
- dh_amplitude: seasonal amplitude in height change [m]
- dh_phase: day of seasonal maximum in height change [day of year]
- dm_trend: trend in mass change [m.w.e./yr]
- dm_amplitude: seasonal amplitude in mass change [m.w.e.]
- dm_phase: day of seasonal maximum in mass change [day of year]
- dm_synthesis_trend: trend in mass change from synthesis of satellite data corrected for changes in fac [m.w.e./yr]
- dm_synthesis_amplitude: seasonal amplitude in mass change from synthesis of satellite data corrected for changes in fac [m.w.e.]
- dm_synthesis_phase: day of seasonal maximum in synthesis of satellite data corrected for changes in fac [day of year]
############################################################################################


############################################################################################
#                         Gardner2025_regional_timseries.nc                                #
############################################################################################
Dimensions:
- rgi: RGI region id (RGI6) 
- date: Center date of monthly data
- error: is 95% confidence interval [false/true]

Variables:
- fac: firn air content anomaly [km^3]
- dm_grace: mass anomaly GRACE [Gt]
- refreeze: refreeze anomaly [Gt]
- rain: rain anomaly [Gt]
- acc: accumulation anomaly [Gt]
- melt: melt anomaly [Gt]
- runoff: runoff anomaly [Gt]
- dm_altim: mass anomaly from firn corrected synthesis of satellite data [Gt]
- dm_glambie: mass anomaly reported by GlaMBIE Team, 2025 [Gt]
- dm: mass anomaly from optimal model fit [Gt]
- dv_altim: volume anomaly from synthesis of satellite data [km^3]
- smb: surface mass balance anomaly from optimal model fit [Gt]
- dv: volume anomaly from optimal model fit [Gt] [km^3]
- ec: evaporation/condensation/deposition/condensation/sublimation anomaly from optimal model fit [Gt]
############################################################################################


############################################################################################
#                          Gardner2025_endorheic_fraction.csv                              #
############################################################################################
Variables:
- rgi: RGI region id (RGI6) 
- region_name: RGI region name (RGI6)        
- dm_endorheic_fraction: fraction of mass change that terminates in a closed basin (i.e. endorheic)
- runoff_endorheic_fraction: fraction of runoff that terminates in a closed basin (i.e. endorheic)
############################################################################################


############################################################################################
#                          Gardner2025_regional_results.csv                                #
############################################################################################
Variables:
- rgi: RGI region id (RGI6) 
- region_name: RGI region name (RGI6)       
- fac_trend_[km^3/yr]: trend in firn air content anomaly ± 95% confidence interval [km^3/yr]
- fac_acceleration_[km^3/yr^2]: acceleration in firn air content anomaly ± 95% confidence interval [km^3/yr^2]
- fac_amplitude_[km^3]: seasonal amplitude in firn air content anomaly ± 95% confidence interval [km^3]
- fac_phase_[day_of_max]: day of maximum in seasonal firn air content anomaly ± 95% confidence interval [day of year]
- refreeze_trend_[Gt/yr]: trend in refreeze anomaly ± 95% confidence interval [Gt/yr]
- refreeze_acceleration_[Gt/yr^2]: acceleration in refreeze anomaly ± 95% confidence interval [Gt/yr^2]
- refreeze_amplitude_[Gt]: seasonal amplitude in refreeze anomaly ± 95% confidence interval [Gt]
- refreeze_phase_[day_of_max]: day of maximum in refreeze anomaly ± 95% confidence interval [day of year]
- rain_trend_[Gt/yr]: trend in rain anomaly ± 95% confidence interval [Gt/yr]
- rain_acceleration_[Gt/yr^2]: acceleration in rain anomaly ± 95% confidence interval [Gt/yr^2]
- rain_amplitude_[Gt]: seasonal amplitude in rain anomaly ± 95% confidence interval [Gt]
- rain_phase_[day_of_max]: day of maximum in seasonal rain anomaly ± 95% confidence interval [day of year]
- acc_trend_[Gt/yr]: trend in accumulation anomaly ± 95% confidence interval [Gt/yr]
- acc_acceleration_[Gt/yr^2]: acceleration in accumulation anomaly ± 95% confidence interval [Gt/yr^2]
- acc_amplitude_[Gt]: seasonal amplitude in accumulation anomaly ± 95% confidence interval [Gt]
- acc_phase_[day_of_max]: day of maximum in accumulation content anomaly ± 95% confidence interval [day of year]
- melt_trend_[Gt/yr]: trend in melt anomaly ± 95% confidence interval [Gt/yr]
- melt_acceleration_[Gt/yr^2]: acceleration in melt anomaly ± 95% confidence interval [Gt/yr^2]
- melt_amplitude_[Gt]: seasonal amplitude in melt anomaly ± 95% confidence interval [Gt]
- melt_phase_[day_of_max]: day of maximum in seasonal melt anomaly ± 95% confidence interval [day of year]
- runoff_trend_[Gt/yr]: trend in runoff anomaly ± 95% confidence interval [Gt/yr]
- runoff_acceleration_[Gt/yr^2]: acceleration in runoff anomaly ± 95% confidence interval [Gt/yr^2]
- runoff_amplitude_[Gt]: seasonal amplitude in runoff anomaly ± 95% confidence interval [Gt]
- runoff_phase_[day_of_max]: day of maximum in seasonal runoff anomaly ± 95% confidence interval [day of year]
- dm_synthesis_trend_[Gt/yr]: trend in mass anomaly from synthesis of satellite data ± 95% confidence interval [Gt/yr]
- dm_synthesis_acceleration_[Gt/yr^2]: acceleration in mass anomaly from synthesis of satellite data ± 95% confidence interval [Gt/yr^2]
- dm_synthesis_amplitude_[Gt]: seasonal amplitude in mass anomaly from synthesis of satellite data ± 95% confidence interval [Gt]
- dm_synthesis_phase_[day_of_max]: day of maximum in seasonal mass anomaly from synthesis of satellite data ± 95% confidence interval [day of year]
- dm_trend_[Gt/yr]: trend in mass anomaly ± 95% confidence interval [Gt/yr]
- dm_acceleration_[Gt/yr^2]: acceleration in mass anomaly ± 95% confidence interval [Gt/yr^2]
- dm_amplitude_[Gt]: seasonal amplitude in mass anomaly ± 95% confidence interval [Gt]
- dm_phase_[day_of_max]: day of maximum in seasonal mass anomaly ± 95% confidence interval [day of year]
- dv_synthesis_trend_[km^3/yr]: trend in volume anomaly from firn corrected synthesis of satellite data ± 95% confidence interval [km^3/yr]
- dv_synthesis_acceleration_[km^3/yr^2]: acceleration in volume anomaly from firn corrected synthesis of satellite data ± 95% confidence interval [km^3/yr^2]
- dv_synthesis_amplitude_[km^3]: seasonal amplitude in volume anomaly from firn corrected synthesis of satellite data ± 95% confidence interval [km^3]
- dv_synthesis_phase_[day_of_max]: day of maximum in seasonal volume anomaly from firn corrected synthesis of satellite data ± 95% confidence interval [day of year]
- smb_trend_[Gt/yr]: trend in surface mass balance anomaly ± 95% confidence interval [Gt/yr]
- smb_acceleration_[Gt/yr^2]: acceleration in surface mass balance anomaly ± 95% confidence interval [Gt/yr^2]
- smb_amplitude_[Gt]: seasonal amplitude in surface mass balance anomaly ± 95% confidence interval [Gt]
- smb_phase_[day_of_max]: day of maximum in seasonal surface mass balance anomaly ± 95% confidence interval [day of year]
- dv_trend_[km^3/yr]: trend in volume anomaly ± 95% confidence interval [km^3/yr]
- dv_acceleration_[km^3/yr^2]: acceleration in volume anomaly ± 95% confidence interval [km^3/yr^2]
- dv_amplitude_[km^3]: seasonal amplitude in volume anomaly ± 95% confidence interval [km^3]
- dv_phase_[day_of_max]: day of maximum in seasonal volume anomaly ± 95% confidence interval [day of year]
- ec_trend_[Gt/yr]: trend in evaporation/condensation/deposition/condensation/sublimation anomaly ± 95% confidence interval [Gt/yr]
- ec_acceleration_[Gt/yr^2]: acceleration in evaporation/condensation/deposition/condensation/sublimation anomaly ± 95% confidence interval [Gt/yr^2]
- ec_amplitude_[Gt]: seasonal amplitude in evaporation/condensation/deposition/condensation/sublimation anomaly ± 95% confidence interval [Gt]
- ec_phase_[day_of_max]: day of maximum in seasonal evaporation/condensation/deposition/condensation/sublimation anomaly ± 95% confidence interval [day of year]
- net_acc_trend_[Gt/yr]: trend in net accumulation anomaly ± 95% confidence interval [Gt/yr]
- net_acc_acceleration_[Gt/yr^2]: acceleration in net accumulation anomaly ± 95% confidence interval [Gt/yr^2]
- net_acc_amplitude_[Gt]: seasonal amplitude in net accumulation anomaly ± 95% confidence interval [Gt]
- net_acc_phase_[day_of_max]: day of maximum in seasonal net accumulation anomaly ± 95% confidence interval [day of year]
- gsi: glacier sustainability index ± 95% confidence interval
- discharge_trend_[Gt/yr]: trend in solid ice discharge through flux gate located near grounding line ± 95% confidence interval [Gt/yr]
- area_km2: glacier area (RGI6) [km^2]
############################################################################################


############################################################################################
#                          Gardner2025_gmax_buffer_population.nc                           #
############################################################################################
Dimensions:
- country: country
- gmax: minimum gmax
- runoff: minimum climatological maximum monthly glacier runoff [m³/s]
- buffer: maximum buffer distance [m]

Variables:
- population: population living within a given buffer distance of a river with a gmax and runoff threshold
############################################################################################


############################################################################################
#                          Gardner2025_glacier_summary_gmax.gpkg                           #
############################################################################################
Variables:
- COMID: MERIT Hydro river reach unique ID
- ocean_terminating: river reach flows into ocean [true/false]
- continent: continent of river reach
- country: continent of river reach
- lengthkm: length of of river reach [km]
- runoff_max_avg: climatological maximum monthly glacier runoff [m³/s]
- runoff_max_month: month of climatological maximum monthly glacier runoff
- gmax_max: maximum gmax between April 1, 2000 and October 1, 2024
- gmax_min: minimum gmax between April 1, 2000 and October 1, 2024
- gmax_range: range in gmax between April 1, 2000 and October 1, 2024
- gmax_avg: average gmax between April 1, 2000 and October 1, 2024
- gmax_month: month of climatological gmax
############################################################################################
