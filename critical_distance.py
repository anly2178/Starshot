import numpy as np

def find_crit_dist(params):
    """
    Input:  The set of parameters
    Output: The critical distance at which the laser begins to spill over the sail. (m)
    """
    #Get parameters
    m_sail = params["m_sail"]
    thickness = params["thickness"]
    density = params["density"]
    k = params["k"]
    laser_size = params["laser_size"]
    wavelength = params["wavelength"]
    alpha = params["alpha"]

    #Calculate the critical distance
    sail_size = np.sqrt(m_sail / (k*density*thickness))
    critical_dist = laser_size * sail_size / (2*wavelength*alpha)
    return critical_dist
