import math
import numpy as np

def string_to_hours(time_string):
    """
    Converts Hours string to float
    :param time_string: Hours String (hh:mm:ss.ss)
    """        
    # Verify separator
    separators = [':', ' ']
    separator = None
    for sep in separators:
        if sep in time_string:
            separator = sep
            break

    if separator is None:
        raise ValueError("Invalid string format. No recognized separator found.")

    components = time_string.split(separator)

    # Check for correct format
    if len(components) != 3:
        raise ValueError(f"Invalid string format. Expected hh{separator}mm{separator}ss.ss")

    hours = abs(int(components[0]))
    minutes = int(components[1])
    seconds = float(components[2])

    total_hours = hours + minutes / 60 + seconds / 3600

    sign = -1 if "-" in time_string else 1
    return sign*total_hours

def string_to_degrees(degrees_string):
    """
    Converts Degrees string to float
    :param degrees_string: Degrees String (dd:mm:ss.ss)
    """
    # Verify separator
    separators = [':', ' ']
    separator = None
    for sep in separators:
        if sep in degrees_string:
            separator = sep
            break

    if separator is None:
        raise ValueError("Invalid string format. No recognized separator found.")

    components = degrees_string.split(separator)

    # Check for correct format
    if len(components) != 3:
        raise ValueError("Invalid string format. Expected dd:mm:ss.ss")

    degrees_int = abs(int(components[0]))
    minutes = int(components[1])    
    seconds = float(components[2])

    degrees = degrees_int + minutes / 60 + seconds / 3600

    sign = -1 if "-" in degrees_string else 1
    return sign*degrees

def calc_ah(ra, lst):
    """Calculates hour angle from RA and Sidereal"""
    ah = 0.2618*(lst - ra)
    if ah > math.pi:
        ah -= 2 * math.pi
    if ah < -math.pi:
        ah += 2 * math.pi
    ah = ah/0.2618
    return ah

def get_az_alt(ra, dec, lst, latitude):
    """Convert equatorial coordinates to horizontal"""
    DEG = 180 / math.pi
    RAD = math.pi / 180.0        
    H = (lst-ra) * 15

    #altitude calc
    sinAltitude = (np.sin(dec * RAD)) * (np.sin(latitude * RAD)) + (np.cos(dec * RAD) * np.cos(latitude * RAD) * np.cos(H * RAD))
    elevation = np.arcsin(sinAltitude) * DEG #altura em graus
    
    #azimuth calc
    y = -1 * np.sin(H * RAD)
    x = (np.tan(dec * RAD) * np.cos(latitude * RAD)) - (np.cos(H * RAD) * np.sin(latitude * RAD))

    #Azimuth calc
    azimuth = np.arctan2(y, x) * DEG
    #converting neg values to pos
    
    if (azimuth.any() < 0):
        azimuth = azimuth + 360 
    if isinstance(azimuth, float):
        if (azimuth < 0):
            azimuth = azimuth + 360
    
    return(azimuth,elevation)

def pol2cart(rho, phi, allsky_angle):
    """Convert polar coordinates to cartesian"""
    """rho: zenith  |   phi: azimuth"""

    # x=f*rho*cos(A-?)+x0
    # y=f*rho*sin(A-?)+y0
    # f = 3.365 piexels/degree
    # rho = star zenith
    # phi = star azimuth
    # allsky_angle = camera rotation
    # x0 e y0 = coordinates of cameras zenith

    x = rho * np.cos(np.radians(phi+allsky_angle))
    y = rho * np.sin(np.radians(phi+allsky_angle))
    if allsky_angle<0:
        x0 = 0.5
        y0 = 12.1
    else:
        x0 = 8.5
        y0 = -5
    
    x=3.365*x+x0
    y=3.365*y+y0
    return(x, y)