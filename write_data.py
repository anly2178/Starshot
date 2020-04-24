from tabulate import tabulate

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
    material = params["material"]
    m_sail = params["m_sail"]
    thickness = params["thickness"]
    area = params["area"]
    density = params["density"]
    reflectivity = params["reflectivity"]
    abs_coeff = params["abs_coeff"]
    absorptance = params["absorptance"]
    k = params["k"]
    power = params["power"]
    laser_size = params["laser_size"]
    wavelength = params["wavelength"]
    alpha = params["alpha"]

    table_params = tabulate([["material","m_sail (kg)","thickness (m)","area (m^2)","density (kgm^-3)","reflectivity","abs_coeff","absorptance","k","power (W)","laser_size (m)","wavelength (m)","alpha"]\
                            ,[material, m_sail,thickness, area, density,reflectivity,abs_coeff,absorptance,k,power,laser_size,wavelength,alpha]])

    #Extract states
    beta = state[0,:]
    dist = state[1,:]
    temp_k = state[2,:]
    # temp_j = state[3,:]
    # table_data = tabulate({"Time (s)": time,"Beta (c)": beta, "Distance (m)": dist, "Kipping Temperature (K)": temp_k, "Justin Temperature (K)": temp_j}, headers="keys", showindex = "always")
    table_data = tabulate({"Time (s)": time,"Beta (c)": beta, "Distance (m)": dist, "Kipping Temperature (K)": temp_k}, headers="keys", showindex = "always")

    f = open(filepath,'w')
    f.write(table_params + '\n\n')
    f.write(table_data)
    f.close()
