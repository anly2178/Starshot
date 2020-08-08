def interpolate_from_list(ls, wavelength):
    """ Fills in any values using a linear fit between data points given in
        the files. Also sets the values beyond the intervals given in the list
        to 0. If the list provided is NoneType or empty, returns 0.
    """

    if ls == None or len(ls) == 0:
        return 0
    
    # From materials, the given list will be in ascending order, so we can
    # use this to our advantage

    if wavelength > ls[-1, 0]:
        return 0        # "worst case scenario in terms of temperature"
    elif wavelength < ls[0, 0]:
        return 0        # also "worst case"
    else:
        # performing linear interpolation between points in the list using the
        # point-gradient formula
        for i, (wl, val) in enumerate(ls):
            print(wl, val)
            if wavelength == wl:
                value = val
                # the first instance where the wavelength input is larger than a
                # value in the list, we want to stop
            elif wavelength < wl:
                # when this is true, we know val lies between the i'th and i-1'th entries
                m = (ls[i, 1] - ls[i-1, 1])/(ls[i, 0] - ls[i-1, 0]) #gradient
                value = m*(wavelength - ls[i-1, 0]) + ls[i-1, 1]
                return value
