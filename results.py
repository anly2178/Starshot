from tabulate import tabulate
import numpy as np

def find_accel_dist_duration(params,state,time):
    """
    From results, finds the acceleration distance (m) and duration (s).
    """
    target = params["target"] #Target speed
    i = 0
    betas = state[0,:]
    dist = state[1,:]
    accel_dist = 0
    while i < len(betas):
        beta = betas[i]
        if beta >= target:
            accel_dist = dist[i]
            accel_dur = time[i]
            return accel_dist, accel_dur
        i += 1

def find_energy_to_launch(params, state, time):
    """
    Finds the energy required to reach the target velocity.
    """
    target = params["target"]
    betas = state[0,:]
    t = 0 #Time taken to reach target velocity
    i = 0
    while i < len(betas):
        beta = betas[i]
        if beta > target:
            t = time[i]
            break
        else:
            i += 1
    power = params["power"]
    energy = power*t
    return energy

def write_results(params, state, time, filepath):
    """
    Writes data to a txt file in tabulated form.
    The data is given in the form of a state and time, where:
        state[0,:] = speed
        state[1,:] = distance
        state[2,:] = temperatures
    Data is formatted as:
        Time (s) | Beta (c) | Distance (m) | Temperature (K)
    """
    #For formatting, paramters are separated into three sections: Sail, optical properties, laser, conditions
    keys = [[],[],[],[]]
    values = [[],[],[],[]]
    #Create index to track which one to add to
    index = 0
    for key, value in params.items():
        #Insert the units
        if key == 'm_sail':
            key += ' (kg)'
        elif key == 'thickness':
            key += ' (m)'
        elif key == 'area':
            key += ' (m^2)'
        elif key == 'radius':
            key += ' (m)'
        elif key == 'density':
            key += ' (kgm^-3)'
        elif key == 'abs_coeff':
            key += ' (cm^-1)'
            index += 1 #Move to optical properties section of parameters
        elif key == 'power':
            key += ' (W)'
            index += 1 #Move to laser section of parameters
        elif key == 'diameter':
            key += ' (m)'
        elif key == 'wavelength':
            key += ' (m)'
        elif key == 'target':
            index += 1 #Move to conditions section of parameters
        elif key == 'accel_dist':
            key += ' (m)'
        elif key == 'max_temp':
            key += ' (K)'
        keys[index].append(key)
        values[index].append(value)
    sail = [keys[0],values[0]]
    optic = [keys[1],values[1]]
    laser = [keys[2],values[2]]
    cond = [keys[3],values[3]]
    table_sail = tabulate(sail)
    table_optic = tabulate(optic)
    table_laser = tabulate(laser)
    table_cond = tabulate(cond)

    #Tabulate values of the state
    beta = state[0,:]
    dist = state[1,:]

    if len(state) == 2:
        table_data = tabulate({"Time (s)": time,"Beta (c)": beta, "Distance (m)": dist}, headers="keys", showindex = "always")
    elif len(state) == 3:
        temp = state[2,:]
        table_data = tabulate({"Time (s)": time,"Beta (c)": beta, "Distance (m)": dist, "Temperature (K)": temp}, headers="keys", showindex = "always")
    #Write results into file
    f = open(filepath,'w')
    f.write(table_sail + '\n' + table_optic + '\n' + table_laser + '\n' + table_cond + '\n\n')
    f.write(table_data)
    f.close()

def strlist_to_float(strlist):
    """
    Converts list of strings to float, for elements that are numeric.
    If it is not a numerical value, keep it as a string.
    """
    i = 0
    while i < len(strlist):
        word = strlist[i]
        try:
            strlist[i] = float(word)
        except ValueError:
            strlist[i] = word
        i += 1
    return strlist

def remove_brackets(ls):
    """
    Given a list, returns a new list with the elements that include brackets
     (  ) removed.
    """
    new_ls = []
    i = 0
    while i < len(ls):
        elem = ls[i]
        if not '(' in elem:
            new_ls.append(elem)
        i += 1
    return new_ls

def trim_2d_zeros(two_array):
    """
    Trim zero columns from 2d array.
    """
    a = two_array
    idx = np.argwhere(np.all(a[..., :] == 0, axis=0))
    a2 = np.delete(a, idx, axis=1)
    return a2

def read_results(filepath):
    """
    Reads data from a txt file, and returns a dictionary of parameters,
    array of times and corresponding multidimensional array of the state.
        state[0,:] = speed
        state[1,:] = distance
        state[2,:] = temperatures
    Might be a bit broken because of the number of parameters causing the
    values to overflow the line. Change to adapt to this in the future.
    """
    #Open file for reading
    f = open(filepath, 'r')

    #Set up lists to insert results
    keys = []
    values = []
    time = np.zeros(200)
    beta = np.zeros(200)
    distance = np.zeros(200)
    temperature = np.zeros(200)
    has_temp = False

    #Introduce index to track what has already been extracted
    index = 0

    while True:
        line = f.readline()
        #Break at end of text file
        if line == '':
            break
        #Skip lines with '---' borders and empty lines
        elements = line.split()
        if len(elements) == 0:
            continue
        elif '-' in elements[0]:
            continue
        #Extracts values according to index
        if index == 0 or index == 2 or index == 4 or index == 6:
            keys += elements
            index += 1
        elif index == 1 or index == 3 or index == 5 or index == 7:
            values += strlist_to_float(elements)
            index += 1
        elif index == 8:
            index += 1 #Reaches state headings
        else:
            #elements is now list of numerical values corresponding to the state
            elements = strlist_to_float(elements)
            i = int(elements[0]) #This is the index for the state arrays
            time[i] = elements[1]
            beta[i] = elements[2]
            distance[i] = elements[3]
            if len(elements) == 5:
                temperature[i] = elements[4]
                has_temp = True
    f.close()
    #Create dictionary of parameters
    keys = remove_brackets(keys)
    params = {}
    i = 0
    while i < len(keys):
        key = keys[i]
        value = values[i]
        params[key] = value
        i += 1
    #Create state multidimensional array
    if not has_temp:
        state = np.array([beta, distance])
    elif has_temp:
        state = np.array([beta,distance,temperature])
    state[:,0] = 1 #To bypass trimming
    state = trim_2d_zeros(state)
    state[:,0] = 0 #Return to zero
    time[0] = 1 #To bypass trimming
    time = np.trim_zeros(time)
    time[0] = 0

    return params, state, time
