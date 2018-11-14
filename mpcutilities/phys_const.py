# /mpcutilities/mpcutilities/phys_const.py
"""
     
 --------------------------------------------------------------
 
 Oct 2018
 
 Payne

 Functions related to
 
 (i) Physical Constants
 
 (ii) Unit-Conversions

 If efficient, might want to replace some of this with the JPL/SPICE stuff
 
 Or novas
 
 Or astropy
 
 Or whatever, doesn't matter ...

 --------------------------------------------------------------
 
"""

# Import third-party packages
# --------------------------------------------------------------
import numpy as np


# Define assorted useful physical constants
# - Might be nicer / tidier to have them in a class / dict
# --------------------------------------------------------------

# Gaussian gravitational constant^2
# 0.01720209895^2 == 2.959122e-04
GMsun          = 2.9591220828559115e-04

# This is now a definition
au_km          = 149597870.700                 
speed_of_light = 2.99792458e5 * 86400./au_km

# Tolerance for light-travel-time calculations
# 1e-2 seconds -> 1e-7 days
LLT_Tol        = 1e-2 / (3600*24)


# Obliquity of ecliptic at J2000
ecl            = (84381.4118*(1./3600)*np.pi/180.)

def rotate_matrix(ecl):
    '''
    Rotation matrices (general) between ecliptic and equatorial coordinates
    '''
    ce = np.cos(ecl)
    se = np.sin(-ecl)
    rotmat = np.array([[1.0, 0.0, 0.0],
                       [0.0,  ce,  se],
                       [0.0, -se,  ce]])
    return rotmat

# Rotation matrix that will turn input ECLIPTIC coordinates into EQUATORIAL coordinates
rot_mat_ec_to_eq  = rotate_matrix(ecl)

def rotate_ec_to_eq( arr ):
    '''
    Simple function to rotate input ECLIPTIC coordinates into EQUATORIAL coordinates
    '''
    return np.dot(rot_mat_ec_to_eq , arr )

# Rotation matrix that will turn input EQUATORIAL coordinates into ECLIPTIC coordinates
rot_mat_eq_to_ec  = rotate_matrix(-ecl)

def rotate_eq_to_ec( arr ):
    '''
    Simple function to rotate input EQUATORIAL coordinates into ECLIPTIC coordinates
    '''
    return np.dot(rot_mat_ec_to_eq , arr )


# Radius of the Earth
Rearth_km      = 6378.1363
Rearth_AU      = Rearth_km/au_km

# Hours-to-degrees
hr2deg         = 15.
