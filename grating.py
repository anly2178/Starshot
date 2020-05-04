import numpy as np
from numpy import arcsin

def find_diffracted_angle(order,wavelength,period):
    """
    Given the grating period and wavlength of light, return the angle (in degrees)
    from the normal for a specific order.
    """
    angle = arcsin(order * wavelength / period)
    angle = np.degrees(angle)
    return angle

def find_rt_coefficients(filepath, incident):
    """
    Input:  filepath to txt file (str)
            incident angle in degrees (int or float)
    Output: Multidimensional array of efficiencies with the form:
                [[r-1, r0, r1, rtot],[t-1, t0, t1, ttot]]
            where r-1 is the reflection coefficient for the -1 order,
            t0 is the transmission coefficient for the 0 order,
            etc.
            and rtot and ttot are the total coefficients.

    Currently specific to the txt file created by Evan Xie, which are the
    results of 2019 Ilic.
    """
    f = open(filepath, 'r')
    #Create structures to store coefficients
    #In order r0,r-1,r1,rtot,t0,t-1,t1,ttot
    coeffs = np.zeros(8)
    #txt file results given in the order r0,r-1,r1,rtot,t0,t1,t-1,ttot
    #Track this with index
    i = 0
    while True:
        line = f.readline()
        if line == '':
            break
        #Split the line into characters
        chars = line.split()
        #Check if the first character is numerical
        try:
            angle = float(chars[0])
        except ValueError:
            continue    #Continue to next line if not numerical
        #If it is numerical, find the desied incident angle
        #The argument incident angle must exactly match the one in the txt file.
        #This can be changed to the closest matching incident angle, but
        #it is unnecessary for now, as we care about normal incidence.
        if angle == incident:
            coeff = float(chars[1])
        else:
            continue    #Continue to next line if not a match
        coeffs[i] = coeff
        i += 1
        if i == 8:
            break
    f.close()
    #Restructure coefficients
    reflec = [coeffs[1],coeffs[0],coeffs[2],coeffs[3]]
    transm = [coeffs[6],coeffs[4],coeffs[5],coeffs[7]]
    coeffs = [reflec,transm]
    return coeffs
