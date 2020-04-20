from tabulate import tabulate
from .motion_with_diffraction import *

def write_data(params, filepath):
    """
    Writes data to a txt file in tabulated form.
    Data is formatted as:
        Time (s) | Beta (c) | Distance (m) | Temperature (K)

    Currently there's no temperature, will update soon.
    """
    #Extract parameters
    m_sail = params["m_sail"]
    thickness = params["thickness"]
    density = params["density"]
    reflectivity = params["reflectivity"]
    k = params["k"]
    power = params["power"]
    laser_size = params["laser_size"]
    wavelength = params["wavelength"]
    alpha = params["alpha"]

    table_params = tabulate([["m_sail (kg)","thickness (m)","density (kgm^-3)","reflectivity","k","power (W)","laser_size (m)","wavelength (m)","alpha"]\
                            ,[m_sail,thickness,density,reflectivity,k,power,laser_size,wavelength,alpha]])

    #Solve for state vs time
    x,t = with_diff_state_vs_t(params)
    beta = x[0,:]
    dist = x[1,:]

    table_data = tabulate({"Time (s)": t,"Beta (c)": beta, "Distance (m)": dist}, headers="keys", showindex = "always")
    f = open(filepath,'w')
    f.write(table_params + '\n\n')
    f.write(table_data)
    f.close()
