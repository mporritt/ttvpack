# ===WIP===
# Functions to fit and optimise for a model including TTVs.

import numpy as np
import lightkurve as lk
import exoplanet as xo
import pymc3 as pm
import pymc3_ext as pmx
from astropy.timeseries import BoxLeastSquares
from astropy.time import Time

from retrieval import get_target
from TimingModel import TimingModel

###
import matplotlib.pyplot as plt
###


def fit_target_tts(lightcurves, **kwargs):  ### MAIN FUNCTION THAT YOU RUN AND IT GIVES YOU TTVS
    """ take a lighkurve collection and fit ttvs
    """
    
    # Group any light curves from consecutive campaigns
    grouped_lcs = [[ lightcurves[0] ]]
    
    for lc in lightcurves[1:]:
        if lc.mission == grouped_lcs[-1][-1].mission:
            sector = 'SECTOR' if (lc.mission == 'TESS') else 'CAMPAIGN'

            if lc.meta[sector] == 1 + grouped_lcs[-1][-1].meta[sector]:
                grouped_lcs[-1].append(lc)
            else:
                grouped_lcs.append([lc])
        else:
            grouped_lcs.append([lc])
            
    # Stitch the grouped light curves
    final_light_curves= [group[0].append(group[1:]) for group in grouped_lcs]
    collection = lk.LightCurveCollection(final_light_curves)
    
    if len(collection) < len(lightcurves):
        print("Target was observed in consecutive sectors/campaigns, stitching adjacent light curves.")
        
    
    # Run fit_ttvs for each grouped light curve
    timing_models = []
    for i, lc in enumerate(collection):
        print(f"Fitting a timing model for data set {i}...")
        timing_models.append( fit_transit_times(lc, **kwargs) )
        
    combined_model = timing_models[0].stitch(timing_models[1:])
    
    return combined_model

    
    
def fit_transit_times(lightcurve, period=None, num_planets=1, MCMC=False):
    """ take a lightcurve and fit ttvs
    """
    if period is not None:
        num_planets = max(num_planets, len(period))
    
    texp = lightcurve.meta['TIMEDEL']    # Exposure time of data
    
    # Specify the flattening window  ### Change to a different method?
    win_len = int( 0.85/texp ) # window of length 0.85 days
    win_len_odd = win_len - 1 + win_len%2
    
    # Flatten and vet the data
    lightcurve = lightcurve\
        .remove_nans()\
        .flatten(window_length=win_len_odd)\
        .remove_outliers(sigma_lower=20, sigma_upper=5) 
        ### make sure this is ok for all data!
    
    t0, period, expected_transit_times = find_transits(lightcurve, period=period, num_planets=num_planets)
    
    t = lightcurve.time.value            # time array
    y = lightcurve.flux.value            # flux array
    
    with pm.Model() as model:
        
        mean = pm.Normal("mean", mu=1.0, sd=0.1) # mean flux
        u = xo.distributions.QuadLimbDark("u") # limb-darkening
        
        r = pm.Uniform("r", lower=0.01, upper=0.1, testval=0.05, shape=num_planets)
        b = xo.distributions.ImpactParameter("b", ror=r, testval=0.5, shape=num_planets)
        
        # Make a parameter of transit times for each planet
        transit_times = []
        for i in range(num_planets):
            transit_times.append(
                pm.Normal(
                    f"tts_{i}",
                    mu=expected_transit_times[i],
                    sd=1.0,
                    shape=len(expected_transit_times[i]),
                )
            )
        
        orbit = xo.orbits.TTVOrbit(b=b, transit_times=transit_times)
        
        light_curves = xo.LimbDarkLightCurve(u).get_light_curve(
            orbit=orbit, r=r, t=t, texp=texp
        )
        #pm.Deterministic("light_curves", light_curves)
        
        net_light_curve = pm.math.sum(light_curves, axis=-1) + mean
        pm.Normal("obs", mu=net_light_curve, sd=5e-4, observed=y) ### <--- CHANGE sd TO THE FLUX_ERR
        
        map_soln = model.test_point
        map_soln = pmx.optimize(start=map_soln, vars=transit_times)
        map_soln = pmx.optimize(start=map_soln, vars=[r, b])
        map_soln = pmx.optimize(start=map_soln, vars=transit_times)
        map_soln = pmx.optimize(start=map_soln)
        
        ### optional MCMC
        
        tts = [Time(map_soln['tts_'+str(i)], format=lightcurve.time.format) 
               for i in range(num_planets)
              ]
        
    return TimingModel(tts=tts, 
                       lightcurve=lightcurve,
                       **{k:map_soln[k] for k in ('mean', 'u', 'r', 'b')}, 
                      )
    
    
    
def find_transits(lightcurve, period=None, num_planets=1, frequency_factor=None):
    """ Find transit signals in a light curve.
    
    Parameters
    ----------
    lightcurve : lightkurve.LightCurve object
    period : list of float
        A list of periods to look for transit signals around.
    num_planets : int
        The number of planetary signals to look for. If periods are given, will instead look for that many.
        
    Returns
    ----------
    A tuple of (t0, period, transit_times), each parameter being a list of length num_planets.
    """
    
    if num_planets < 1:
        return ([], [], [])    # Recursion end
    
    if period is not None:
        num_planets = max(num_planets, len(period))
    else:
        period = []
        
    # Run BLS for the first planet
    
    bls = BoxLeastSquares(lightcurve.time.value, lightcurve.flux.value)    
    
    
    if frequency_factor is None:
        ff = 1.0 * (lightcurve.time[-1] - lightcurve.time[0]).value / 75 
    else:
        ff = frequency_factor
    sample_durations = [0.05,0.1,0.2]
    
    if len(period) != 0:
        sample_period = period[0]
        sample_periods = bls.autoperiod(sample_durations, 
                                        minimum_period=sample_period*0.95, 
                                        maximum_period=sample_period*1.05, 
                                        frequency_factor=ff/2)
    else:
        sample_periods = bls.autoperiod(sample_durations, frequency_factor=ff)
        
    periodogram = bls.power(sample_periods, sample_durations, objective="snr")
    
    # ###
    # plt.figure()
    # plt.plot(periodogram.period, periodogram.power)
    # plt.title("Periodogram")
    # plt.show()
    # ###
        
    max_power  = np.argmax(periodogram.power)
    new_t0     = periodogram.transit_time[max_power]
    new_period = periodogram.period[max_power]
    duration   = periodogram.duration[max_power]
    
    new_transit_times = xo.orbits.ttv.compute_expected_transit_times(
        lightcurve.time.value[0], lightcurve.time.value[-1], new_period, new_t0)[0]
    
    ### Use bls..compute_stats() to find out more and vet planets???
    
    print(f"BLS found signal with t0={new_t0:.2f}, period={new_period:.4f}, num_transits={len(new_transit_times)}")
    
    # Mask the found transits
    mask = lightcurve.create_transit_mask(new_period, new_t0, duration*2)
    
    # Recursive call for additional planets
    rec_t0, rec_period, rec_transit_times = find_transits(lightcurve[np.invert(mask)], 
                                              period=period[1:], 
                                              num_planets=num_planets-1
                                             )
    
    return ([new_t0] + rec_t0, [new_period] + rec_period, [new_transit_times] + rec_transit_times)