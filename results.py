from tabulate import tabulate
import numpy as np

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
    #Tabulates the parameters
    keys = []
    values = []
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
        elif key == 'power':
            key += ' (W)'
        elif key == 'diameter':
            key += ' (m)'
        elif key == 'wavelength':
            key += ' (m)'
        keys.append(key)
        values.append(value)
    table = [keys,values]
    table_params = tabulate(table)

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
    f.write(table_params + '\n\n')
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

def read_results(filepath):
    """
    Reads data from a txt file, and returns a dictionary of parameters,
    array of times and corresponding multidimensional array of the state.
        state[0,:] = speed
        state[1,:] = distance
        state[2,:] = temperatures
    """
    #Open file for reading
    f = open(filepath, 'r')

    #Set up lists to insert results
    keys = []
    values = []
    states = []
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
        if index == 0:
            keys = elements
            index += 1
        elif index == 1:
            values = strlist_to_float(elements)
            index += 1
        elif index == 2:
            states = elements
            index += 1
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

    return params, state, time
