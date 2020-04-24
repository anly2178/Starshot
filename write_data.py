from tabulate import tabulate

def write_data(params, state, time, filepath):
    """
    Writes data to a txt file in tabulated form.
    The data is given in the form of a state and time, where:
        state[0,:] = speed
        state[1,:] = distance
        state[2,:] = temperatures
    Data is formatted as:
        Time (s) | Beta (c) | Distance (m) | Temperature (K)
    """
    #Extract parameters
    material = params["material"]
    m_sail = params["m_sail"]
    thickness = params["thickness"]
    area = params["area"]
    density = params["density"]
    reflectance = params["reflectance"]
    abs_coeff = params["abs_coeff"]
    absorptance = params["absorptance"]
    power = params["power"]
    diameter = params["diameter"]
    wavelength = params["wavelength"]

    table_params = tabulate([["material","m_sail (kg)","thickness (m)","area (m^2)","density (kgm^-3)","reflectance","abs_coeff (cm^-1)","absorptance","power (W)","diameter (m)","wavelength (m)"]\
                            ,[material, m_sail,thickness, area, density,reflectance,abs_coeff,absorptance,power,diameter,wavelength]])

    #Extract states
    beta = state[0,:]
    dist = state[1,:]
    temp = state[2,:]

    table_data = tabulate({"Time (s)": time,"Beta (c)": beta, "Distance (m)": dist, "Temperature (K)": temp}, headers="keys", showindex = "always")

    f = open(filepath,'w')
    f.write(table_params + '\n\n')
    f.write(table_data)
    f.close()
