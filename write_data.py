from tabulate import tabulate
from .motion_with_diffraction import *

def write_data(params, state, time, filepath):
    """
    Writes data to a txt file in tabulated form.
    The data is given in the form of a state and time, where:
        state[0,:] = speed
        state[1,:] = distance
        state[2,:] = temperatures (Optional)
    Data is formatted as:
        Time (s) | Beta (c) | Distance (m) | Temperature (K)
    """
    #Extract parameters
    m_sail = params["m_sail"]
    thickness = params["thickness"]
    density = params["density"]
    reflectivity = params["reflectivity"]
    absorptance = params["absorptance"]
    k = params["k"]
    power = params["power"]
    laser_size = params["laser_size"]
    wavelength = params["wavelength"]
    alpha = params["alpha"]

    table_params = tabulate([["m_sail (kg)","thickness (m)","density (kgm^-3)","reflectivity","absorptance","k","power (W)","laser_size (m)","wavelength (m)","alpha"]\
                            ,[m_sail,thickness,density,reflectivity,absorptance,k,power,laser_size,wavelength,alpha]])

    #Extract states
    beta = state[0,:]
    dist = state[1,:]
    temp = state[2,:]
    table_data = tabulate({"Time (s)": time,"Beta (c)": beta, "Distance (m)": dist, "Temperature (K)": temp}, headers="keys", showindex = "always")

    f = open(filepath,'w')
    f.write(table_params + '\n\n')
    f.write(table_data)
    f.close()
