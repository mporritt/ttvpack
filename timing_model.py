# ===WIP===


class TimingModel(system_parameters, times, times_err):
    """ Stores and presents the results of model optimisation.
    """
    def __init__(self, system_parameters, times, time_err):
        self.system_parameters = system_parameters
        self.times = times
        self.times_err = times_err
        
        
    def visualise(self):
        """ Present several relevent plots.
        """
        
    