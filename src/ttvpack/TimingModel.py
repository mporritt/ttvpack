# ===WIP===

import numpy as np
import matplotlib.pyplot as plt
import itertools
from astropy.timeseries import BoxLeastSquares


__all__ = ["TimingModel"]


class TimingModel(dict):
    """ Stores and presents the optimised ttv model.
    """
    
    def __init__(self, mean=None, u=None, r=None, b=None, tts=None, lightcurve=None):
        
        self.mean = mean
        self.u = u
        self.r = r
        self.b = b
        self.tts = tts
        self.lightcurve = lightcurve
        
        self.is_stitched = False
        self.section_names  = None
        
        
    def ttvs(self, **kwargs):
        if self.is_stitched:
            return ttvs_from_tts(self.tts, **kwargs)
        else:
            return ttvs_from_tts([self.tts], **kwargs)
    
    
    def stitch(self, others):
        ### make sure to stitch the right planets!!!
        
        models = [self] + others
        
        # Copy the parameters into lists with an entry for each sector
        mean = [model.mean for model in models]
        u    = [model.u    for model in models]
        r    = [model.r    for model in models]
        b    = [model.b    for model in models]
        tts  = [model.tts.copy() for model in models]  # List of transit times indexed as tts[sector][planet][transit_num]
        
        # Put the times in the same units
        for sector in range(len(tts)):
            for planet in range(len(tts[sector])):
                tts[sector][planet].format = self.tts[0].format
        
        # Stitch the light curves into one
        lightcurve = self.lightcurve.append([other.lightcurve for other in others])

        
        new_timing_model = self.__class__(mean=mean, u=u, r=r, b=b, tts=tts, lightcurve=lightcurve)
        new_timing_model.is_stitched = True
        
        # Make a list of name strings with the mission and sector/campaign
        new_timing_model.section_names = []
        for model in models:
            if model.lightcurve.mission.lower() == 'k2':
                name = "k2 campaign " + str(model.lightcurve.campaign).zfill(2)
            elif model.lightcurve.mission.lower() == 'tess':
                name = "TESS sector " + str(model.lightcurve.sector).zfill(2)
            else:
                name = "unknown campaign/sector"
            new_timing_model.section_names.append(name)
        
        return new_timing_model
    
    
    def plot_ttvs(self):
        """ Plots an O-C diagram of the transit timing variations.
        """
        ttvs, period, t0 = self.ttvs(return_period=True)
        
        
        plt.figure()
        for planet in range(len(ttvs)):
            planet_tts = [time.value for tts_sector in self.tts for time in tts_sector[planet]]
            plt.plot(planet_tts, 1440*np.array(ttvs[planet]), 'o', label=f"period={period[planet]:.2f}")
        
        plt.axhline(color='k', linestyle=':')
        plt.title("O - C Diagram")
        plt.xlabel(f"Time [{self.tts[0][0].format}]")
        plt.ylabel("Timing Variations [mins]")
        plt.legend()
        plt.show()
        
        
        
    def visualise(self):
        """ Present several relevent plots.
        """
        return
    
    
    def __repr__(self):
        
        if self.is_stitched:
            mean_str = ", ".join([ "{:.3f}".format(mean) for mean in self.mean ])
            u_str    = ", ".join([ "({:.3f}, {:.3f})".format(u[0], u[1]) for u in self.u ])
            
            r_str = [] # list of strings for each planet
            b_str = []
            tts_str = []
            
            for planet in range(len(self.r[0])):
                r_str_planet = [ "{:.3f}".format(r) for r in np.array(self.r)[:,planet] ]
                r_str.append(   "(" + ", ".join(r_str_planet) + ")" )
                
                b_str_planet = [ "{:.3f}".format(b) for b in np.array(self.b)[:,planet] ]
                b_str.append(   "(" + ", ".join(b_str_planet) + ")" )
                
                tts_str_planet = [ "{:.3f}".format(tt.value) for tts_sector in self.tts for tt in tts_sector[planet]]
                tts_str.append( "[" + ", ".join(tts_str_planet) + "]" )
                
        else:
            mean_str = "{:.3f}".format(self.mean)
            u_str = "({:.3f}, {:.3f})".format(self.u[0], self.u[1])
            r_str = [ "{:.3f}".format(r) for r in self.r ]
            b_str = [ "{:.3f}".format(b) for b in self.b ]
            tts_str = [ "[ " + ", ".join([ "{:.3f}".format(t) for t in tts ]) + " ]"
                       for tts in self.tts 
                      ]
            
        
        out = ("\n"
               "========== Model Parameters ========== \n"
              f"mean = {mean_str} \n"
              f"u = {u_str} \n")
        if self.is_stitched:
            out+=(
              f"missions = {self.section_names} \n")
        out += "-------------------------------------- \n"
        
        for i in range(len(self.r[0])):
            out += (
              f"          Planet {i} Parameters          \n"
              f" r = {r_str[i]} r/R* \n"
              f" b = {b_str[i]} \n"
              f" tts = {tts_str[i]} {self.tts[0][0].format} \n"
               "-------------------------------------- \n"
            )
            
        return out      
        
        
        
        
def ttvs_from_tts(transit_times, return_period=False):
    """ Calculate the variations in transit timings.
    Performs a least-squares fit on the data to obtain the reference period and t0.
        
    Parameters
    ----------
    tts : list
        Transit times indexed as transit_times[sector][planet][transit_num].
    return_period : bool, opt.
        If True, returns a tuple of (ttvs, period, t0).
        
    Returns
    ----------
    ttvs : list
        List of ttvs, with an entry for each planet.
    """
    
    return_ttvs = []
    return_periods = []
    return_t0s = []
    
    for planet in range(len(transit_times[0])):
        tts = [tts_sec[planet].value for tts_sec in transit_times]
        tts_flat = [val for tts_sec in tts for val in tts_sec]

        # Get approximate period
        period_est = np.median(np.diff(tts_flat))

        # for each gap in data, get an approximate number of transits in between, and then get a list of gap sizes to try for each data gap
        sample_gap_sizes = []
        for i in range(len(tts)-1):
            gap = (tts[i+1][0] - tts[i][-1] - period_est/2) // period_est
            sample_range = np.ceil(0.05*gap + 1)
            sample_array = np.arange(gap-sample_range, gap+sample_range+1, dtype=int)
            sample_gap_sizes.append(sample_array)

        combinations = list(itertools.product(*sample_gap_sizes)) # all combinations of different gap sizes


        best_ttvs = None
        best_transit_numbers = None
        best_t0 = 0
        best_period = 0
        best_error = float('inf')

        for gaps in combinations:

            # Determine the transit numbers of each transit time given the gaps between sectors
            transit_numbers = np.arange(len(tts_flat))
            for i in range(len(gaps)):
                ind = sum([ len(tts_sec) for tts_sec in tts[:i+1] ])
                transit_numbers[ind:] += gaps[i]

            # Fit a line to the data and get the ttvs
            period, t0 = np.polyfit(transit_numbers, tts_flat, 1)

            expected_transit_times = [t0 + n * period for n in transit_numbers]
            ttvs = np.subtract(tts_flat, expected_transit_times)

            sum_squared_err = np.sum(np.square(ttvs))

            if sum_squared_err < best_error:
                best_ttvs = ttvs
                best_transit_numbers = transit_numbers
                best_t0 = t0
                best_period = period
                best_error = sum_squared_err
                
        return_periods.append(best_period)
        return_t0s.append(best_t0)
        return_ttvs.append(best_ttvs)
            
    if return_period:
        return (return_ttvs, return_periods, return_t0s)
    return return_ttvs

#    ttvs = []
#    periods = []
#    t0s = []
#    
#
#    for planet in range(len(transit_times[0])):
#        tts = [val for tts_sec in transit_times for val in tts_sec[planet]]
#        
#        transit_numbers = _get_transit_numbers(tts)
#        
#        period, t0 = np.polyfit(transit_numbers, tts, 1)
#
#        expected_transit_times = [t0 + n * period for n in transit_numbers]
#        ttv = tts - expected_transit_times
#
#        ttvs.append(ttv)
#        periods.append(period)
#        t0s.append(t0s)
#
#    if return_period:
#        return (ttvs, periods, t0s)
#
#    return ttvs


def _get_transit_numbers(tts):
    """ Determines the transit number of each transit given a list of known transit times.
    Should identify gaps in the transit times associated with stitched data sets.
    
    Parameters
    ----------
    tts : list of float
        List of transit times for one planet.
        
    Returns
    ----------
    transit_numbers : list of int
        A list of transit numbers the same length as tts.
    """
    
    period_est = np.median(np.diff(tts))
    
    linear_transit_times = np.arange(tts[0] - period_est, tts[-1] + period_est, period_est)

    best_fit_dt = 0
    best_fit_error = float("inf")

    # Scan the linear times accross to find the offset with the least sum-square-error.
    for dt in np.arange(0, period_est.value, 1/1440):
        print('bop', end='')
        sum_error = 0

        for tt in tts:
            dist = np.abs(np.subtract(tt, linear_transit_times + dt))
            error = np.square( np.min(dist) )
            sum_error += error

        if sum_error < best_fit_error:
            best_fit_dt = dt
            best_fit_error = sum_error

    # now retrieve which entries correspond to which linear transit times
    transit_numbers = []
    for tt in tts:
        dist = np.abs(np.subtract(tt, linear_transit_times + best_fit_dt))
        transit_numbers.append( np.argmin(dist) )
    
    # offset to make the first transit the zeroth
    transit_numbers = [n - transit_numbers[0] for n in transit_numbers]

    return transit_numbers
    