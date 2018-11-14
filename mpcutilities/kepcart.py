"""

 --------------------------------------------------------------

 Oct 2018
 
 Payne

 Derived from Holman's previous kepcart code

 Use C to do fast coordinate conversions (Keplerian <-> Cartesian)

 --------------------------------------------------------------
 
"""

# Import third-party packages
# --------------------------------------------------------------
import os
import numpy as np
from ctypes import *
from pkg_resources import resource_filename

# Importing of local modules/packages
# --------------------------------------------------------------
import mpcutilities.classes     as Classes

# Import local files / dirs  
# --------------------------------------------------------------
#lib = CDLL(os.path.join(os.path.dirname(__file__), 'kepcart_src/libkepcart.so'))
#lib = CDLL(resource_filename('mpcutilities','kepcart.so'))
#lib = CDLL('libkepcart.so')
lib = CDLL(os.path.join(os.path.dirname(__file__), 'libkepcart.so'))


# Define "kepcart" routines
# --------------------------------------------------------------
def cart2kep(GM, cartState):
    """
    Converts cartesian coordinates to keplerian coordinates
    
    Keplerians elements are: (a, e, incl, longnode, argperi, meananom)
    
    ***CONVERSION USES A CALL TO C-CODE***
    
    Parameters
    ----------
    GM          : float,
        Constant
        
    cartState   : "CartState" Object-type as defined in MPCFormat.
        Assumes HELIOCENTRIC ECLIPTIC CARTESIAN initial conditions
    
    Returns
    -------
    (a, e, incl, longnode, argperi, meananom) : tuple of floats
    
    Examples
    --------
    >>> ...
    
    """


    _cart2kep = lib.cart2kep
    _cart2kep.argtypes = (c_double, Classes.CartState)
    _cart2kep.restype = None

    a = c_double()
    e = c_double()
    incl = c_double()
    longnode = c_double()
    argperi = c_double()
    meananom = c_double()

    return_value = _cart2kep(GM, cartState, byref(a), byref(e), byref(incl), byref(longnode), byref(argperi), byref(meananom))

    return (a.value, e.value, incl.value, longnode.value, argperi.value, meananom.value)


def keplerian(GM, cartState):
    """
    Identical to cart2kep
    
    Converts cartesian coordinates to keplerian coordinates
    
    Provided so that Holman's legacy code will always work
    """
    #   return cart2kep(GM, cartState)

    _keplerian = lib.keplerian
    _keplerian.argtypes = (c_double, Classes.CartState)
    _keplerian.restype = None
    
    a = c_double()
    e = c_double()
    incl = c_double()
    longnode = c_double()
    argperi = c_double()
    meananom = c_double()
    
    return_value = _keplerian(GM, cartState, byref(a), byref(e), byref(incl), byref(longnode), byref(argperi), byref(meananom))
    
    return (a.value, e.value, incl.value, longnode.value, argperi.value, meananom.value)






def cart2kep_array(GM, cartStateArray):
    """
    Converts arrays of cartesian coordinates to arrays of keplerian coordinates
    
    Keplerians elements are: (a, e, incl, longnode, argperi, meananom)
    
    ***CONVERSION USES A CALL TO C-CODE***
    
    Parameters
    ----------
    GM          : float,
        Constant
    cartState   : "CartStateArray" Object-type as defined in MPCFormat.
    
        Assumes HELIOCENTRIC ECLIPTIC CARTESIAN initial conditions
        Length = N_s
    
    Returns
    -------
    a, e, incl, longnode, argperi, meananom : numpy arrays
    
    Examples
    --------
    >>> ...
    
    """

    num = len(cartStateArray)
    
    StateArray = Classes.CartState * num

    a_arr = np.zeros((num), dtype=np.double)
    e_arr = np.zeros((num), dtype=np.double)
    incl_arr = np.zeros((num), dtype=np.double)
    longnode_arr = np.zeros((num), dtype=np.double)
    argperi_arr = np.zeros((num), dtype=np.double)
    meananom_arr =np.zeros((num), dtype=np.double)

    _cart2kep_array = lib.cart2kep_array
    _cart2kep_array.argtypes = (c_int, c_double, POINTER(StateArray), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double))
    _cart2kep_array.restype = None


    return_value = _cart2kep_array(num, GM, byref(cartStateArray),
                               a_arr.ctypes.data_as(POINTER(c_double)),
                               e_arr.ctypes.data_as(POINTER(c_double)),
                               incl_arr.ctypes.data_as(POINTER(c_double)),
                               longnode_arr.ctypes.data_as(POINTER(c_double)),
                               argperi_arr.ctypes.data_as(POINTER(c_double)),
                               meananom_arr.ctypes.data_as(POINTER(c_double)))

    return a_arr, e_arr, incl_arr, longnode_arr, argperi_arr, meananom_arr

def keplerians(GM, cartStateArray):
    """
    Identical to cart2kep_array
    
    Converts arrays of cartesian coordinates to arrays of keplerian coordinates
    
    Provided so that Holman's legacy code will always work
    """
    num = len(cartStateArray)
    
    StateArray = Classes.CartState * num
    
    a_arr = np.zeros((num), dtype=np.double)
    e_arr = np.zeros((num), dtype=np.double)
    incl_arr = np.zeros((num), dtype=np.double)
    longnode_arr = np.zeros((num), dtype=np.double)
    argperi_arr = np.zeros((num), dtype=np.double)
    meananom_arr =np.zeros((num), dtype=np.double)
    
    _keplerians = lib.keplerians
    _keplerians.argtypes = (c_int, c_double, POINTER(StateArray), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double))
    _keplerians.restype = None
    
    
    return_value = _keplerians(num, GM, byref(cartStateArray),
                                   a_arr.ctypes.data_as(POINTER(c_double)),
                                   e_arr.ctypes.data_as(POINTER(c_double)),
                                   incl_arr.ctypes.data_as(POINTER(c_double)),
                                   longnode_arr.ctypes.data_as(POINTER(c_double)),
                                   argperi_arr.ctypes.data_as(POINTER(c_double)),
                                   meananom_arr.ctypes.data_as(POINTER(c_double)))
                                   
    return a_arr, e_arr, incl_arr, longnode_arr, argperi_arr, meananom_arr







def kep2cartState(GM, a, e, incl, longnode, argperi, meananom):
    """
    Converts keplerian coordinates to a cartesian state
    
    Keplerians elements are: (a, e, incl, longnode, argperi, meananom)
    
    Assumes HELIOCENTRIC ECLIPTIC KEPLERIAN initial conditions
    
    ***CONVERSION USES A CALL TO C-CODE***
    
    Parameters
    ----------
    GM          :   float
        Gravity
    a           :   float
        Semi-major axis
    e           :   float 
        Eccentricity
    incl        :   float 
        Inclination
    longnode    :   float   
        Longitude of ascending node
    argperi     :   float 
        Argument of pericenter
    meananom    :   float
        Mean anomaly
    
    
    Returns
    -------
    cartState   :   "CartState" Object-type as defined in MPCFormat.
        Assumes HELIOCENTRIC ECLIPTIC CARTESIAN
    
    Examples
    --------
    >>> ...
    
    """
    
    _kep2cartState = lib.kep2cartState
    _kep2cartState.argtypes = (c_double, c_double, c_double, c_double, c_double, c_double, c_double, POINTER(Classes.CartState))
    _kep2cartState.restype = None

    cartState = Classes.CartState(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    return_value = _kep2cartState(GM, a, e, incl, longnode, argperi, meananom, byref(cartState))

    return cartState


def cartesian(GM, a, e, incl, longnode, argperi, meananom):
    """
    Identical to kep2cartState
    
    Converts keplerian coordinates to a cartesian state
    
    Provided so that Holman's legacy code will always work
    """
    _cartesian = lib.cartesian
    _cartesian.argtypes = (c_double, c_double, c_double, c_double, c_double, c_double, c_double, POINTER(Classes.CartState))
    _cartesian.restype = None
    
    cartState = Classes.CartState(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    
    return_value = _cartesian(GM, a, e, incl, longnode, argperi, meananom, byref(cartState))
    
    return cartState







def kep2cartStateArray(GM, a_arr, e_arr, incl_arr, longnode_arr, argperi_arr, meananom_arr):
    """
    Converts arrays of keplerian coordinates to a cartesian state array
    
    Keplerians elements are: (a, e, incl, longnode, argperi, meananom)
    
    Assumes HELIOCENTRIC ECLIPTIC KEPLERIAN initial conditions
    
    ***CONVERSION USES A CALL TO C-CODE***
    
    Parameters
    ----------
    GM          :   float
        Gravity
    a           :   array
        Semi-major axis
    e           :   array
        Eccentricity
    incl        :   array
        Inclination
    longnode    :   array
        Longitude of ascending node
    argperi     :   array
        Argument of pericenter
    meananom    :   array
        Mean anomaly
    
    
    Returns
    -------
    (x, y, z, xd, yd, zd) : tuple of floats    <<-- *** SOME BACK TO THIS: WHAT IS RETURNED ??? -- PRESUME IT'S A STATE -- IF SO REPLACE WITH LINE BELOW ***
    
    cartStateArray      :   "CartState" Object-type as defined in MPCFormat.
        Assumes HELIOCENTRIC ECLIPTIC CARTESIAN


    Examples
    --------
    >>> ...
    
    """
    num = len(a_arr)
    StateArray = Classes.CartState * num
    state_arr  = StateArray()

    _kep2cartStateArray = lib.kep2cartStateArray
    _kep2cartStateArray.argtypes = (c_int, c_double, POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(StateArray))
    _kep2cartStateArray.restype = None


    return_value = _kep2cartStateArray(num, GM,
                               a_arr.ctypes.data_as(POINTER(c_double)),
                               e_arr.ctypes.data_as(POINTER(c_double)),
                               incl_arr.ctypes.data_as(POINTER(c_double)),
                               longnode_arr.ctypes.data_as(POINTER(c_double)),
                               argperi_arr.ctypes.data_as(POINTER(c_double)),
                               meananom_arr.ctypes.data_as(POINTER(c_double)),
                               byref(state_arr))

    return state_arr

def cartesians(GM, a_arr, e_arr, incl_arr, longnode_arr, argperi_arr, meananom_arr):
    """
    Identical to kep2cartStateArray
    
    Converts arrays of keplerian coordinates to a cartesian state array
    
    Provided so that Holman's legacy code will always work
    """
    num = len(a_arr)
    StateArray = Classes.CartState * num
    state_arr  = StateArray()
    
    _cartesians = lib.cartesians
    _cartesians.argtypes = (c_int, c_double, POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(StateArray))
    _cartesians.restype = None
    
    
    return_value = _cartesians(num, GM,
                                       a_arr.ctypes.data_as(POINTER(c_double)),
                                       e_arr.ctypes.data_as(POINTER(c_double)),
                                       incl_arr.ctypes.data_as(POINTER(c_double)),
                                       longnode_arr.ctypes.data_as(POINTER(c_double)),
                                       argperi_arr.ctypes.data_as(POINTER(c_double)),
                                       meananom_arr.ctypes.data_as(POINTER(c_double)),
                                       byref(state_arr))
                                       
    return state_arr




def kep2cartPV(GM, a_arr, e_arr, incl_arr, longnode_arr, argperi_arr, meananom_arr):
    """
    Converts arrays of keplerian coordinates to arrays of cartesian positions and velocities
    
    Keplerians elements are: (a, e, incl, longnode, argperi, meananom)
    
    Assumes HELIOCENTRIC ECLIPTIC KEPLERIAN initial conditions
    
    ***CONVERSION USES A CALL TO C-CODE***
    
    Parameters
    ----------
    GM          :   float
        Gravity
    a           :   array
        Semi-major axis
    e           :   array
        Eccentricity
    incl        :   array
        Inclination
    longnode    :   array
        Longitude of ascending node
    argperi     :   array
        Argument of pericenter
    meananom    :   array
        Mean anomaly
    
    
    Returns
    -------
    pos_arr, vel_arr    : ndarrays
        ** DESCRIBE THE ORDER THAT THESE ARE IN ***
    
    
    Examples
    --------
    >>> ...
    
    """

    num  = len(a_arr)
    size = num*3
    array_of_size_doubles = c_double*size

    pos_arr = array_of_size_doubles()
    vel_arr = array_of_size_doubles()

    _kep2cartPV = lib.kep2cartPV
    _kep2cartPV.argtypes = (c_int, c_double, POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(array_of_size_doubles), POINTER(array_of_size_doubles))
    _kep2cartPV.restype = None


    return_value = _kep2cartPV(num, GM,
                                      a_arr.ctypes.data_as(POINTER(c_double)),
                                      e_arr.ctypes.data_as(POINTER(c_double)),
                                      incl_arr.ctypes.data_as(POINTER(c_double)),
                                      longnode_arr.ctypes.data_as(POINTER(c_double)),
                                      argperi_arr.ctypes.data_as(POINTER(c_double)),
                                      meananom_arr.ctypes.data_as(POINTER(c_double)),
                                      byref(pos_arr),
                                      byref(vel_arr)
                               )
                               
    # At this stage the output is flat:
    # E.g. x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3
    # Seems of general use to reshape it ...
    XYZ = np.array(pos_arr).reshape((-1,3))
    UVW = np.array(vel_arr).reshape((-1,3))
    
    return XYZ, UVW

def cartesian_vectors(GM, a_arr, e_arr, incl_arr, longnode_arr, argperi_arr, meananom_arr):
    """
    Identical to kep2cartPV
    
    Converts arrays of keplerian coordinates to arrays of cartesian positions and velocities
    
    Provided so that Holman's legacy code will always work
    """
    num  = len(a_arr)
    size = num*3
    array_of_size_doubles = c_double*size
    
    pos_arr = array_of_size_doubles()
    vel_arr = array_of_size_doubles()
    
    _cartesian_vectors = lib.cartesian_vectors
    _cartesian_vectors.argtypes = (c_int, c_double, POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(array_of_size_doubles), POINTER(array_of_size_doubles))
    _cartesian_vectors.restype = None
    
    
    return_value = _cartesian_vectors(num, GM,
                               a_arr.ctypes.data_as(POINTER(c_double)),
                               e_arr.ctypes.data_as(POINTER(c_double)),
                               incl_arr.ctypes.data_as(POINTER(c_double)),
                               longnode_arr.ctypes.data_as(POINTER(c_double)),
                               argperi_arr.ctypes.data_as(POINTER(c_double)),
                               meananom_arr.ctypes.data_as(POINTER(c_double)),
                               byref(pos_arr),
                               byref(vel_arr)
                               )
        
    # At this stage the output is flat:
    # E.g. x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3
    # Seems of general use to reshape it ...
    XYZ = np.array(pos_arr).reshape((-1,3))
    UVW = np.array(vel_arr).reshape((-1,3))
    
    return XYZ, UVW






def kepState2cartPV(GM, elementsArray):
    """
    Converts arrays of keplerian element-objects to arrays of cartesian positions and velocities
    Assumes HELIOCENTRIC ECLIPTIC KEPLERIAN initial conditions
    
    ***CONVERSION USES A CALL TO C-CODE***
    
    Parameters
    ----------
    GM              :   float
        Gravity
    elementsArray   :   "elementsArray" Object-type as defined in MPCFormat.
    
    
    Returns
    -------
    pos_arr, vel_arr    :
        ** DESCRIBE THESE ***
    
    
    Examples
    --------
    >>> ...
    
    """
    num             = len(elementsArray)
    ElementsArray   = Classes.KepState * num
    
    size = num*3
    array_of_size_doubles = c_double*size

    pos_arr = array_of_size_doubles()
    vel_arr = array_of_size_doubles()

    _kepState2cartPV = lib.kepState2cartPV
    _kepState2cartPV.argtypes = (c_int, c_double, POINTER(ElementsArray), POINTER(array_of_size_doubles), POINTER(array_of_size_doubles))
    _kepState2cartPV.restype = None

    return_value = _kepState2cartPV(num, GM,
                                    byref(elementsArray),
                                    byref(pos_arr),
                                    byref(vel_arr))
                                    
                                    
    # At this stage the output is flat:
    # E.g. x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3
    # Seems of general use to reshape it ...
    XYZ = np.array(pos_arr).reshape((-1,3))
    UVW = np.array(vel_arr).reshape((-1,3))
    
    return XYZ, UVW

def cartesian_elements(GM, elementsArray):
    """
    Identical to kepState2cartPV
    
    Converts arrays of keplerian element-objects to arrays of cartesian positions and velocities
    
    Provided so that Holman's legacy code will always work
    """
    num             = len(elementsArray)
    ElementsArray   = Classes.KepState * num
    
    size = num*3
    array_of_size_doubles = c_double*size
    
    pos_arr = array_of_size_doubles()
    vel_arr = array_of_size_doubles()
    
    _cartesian_elements = lib.cartesian_elements
    _cartesian_elements.argtypes = (c_int, c_double, POINTER(ElementsArray), POINTER(array_of_size_doubles), POINTER(array_of_size_doubles))
    _cartesian_elements.restype = None
    
    return_value = _cartesian_elements(num, GM,
                                    byref(elementsArray),
                                    byref(pos_arr),
                                    byref(vel_arr))
        
    # At this stage the output is flat:
    # E.g. x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3
    # Seems of general use to reshape it ...
    XYZ = np.array(pos_arr).reshape((-1,3))
    UVW = np.array(vel_arr).reshape((-1,3))
    
    return XYZ, UVW



