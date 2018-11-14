"""

--------------------------------------------------------------

Nov 2018
 - Payne

Holman's WIP to generalize kepcart
 - *** INCOMPLETE & UNTESTED ***

Use C to do fast coordinate conversions (Keplerian <-> Cartesian)


--------------------------------------------------------------

"""

from ctypes import *
import numpy as np
from numpy.ctypeslib import ndpointer

class State(Structure):
    _fields_ = [
        ('x', c_double),
        ('y', c_double),
        ('z', c_double),
        ('xd', c_double),
	('yd', c_double),
	('zd', c_double)
    ]

class Elements(Structure):
    _fields_ = [
        ('a', c_double),
        ('e', c_double),
        ('incl', c_double),
        ('longnode', c_double),
	('argperi', c_double),
	('meananom', c_double)
    ]

#lib = CDLL('/Users/mholman/Dropbox/MPC/libkepcart.so')
#lib = CDLL("/Users/matthewjohnpayne/Dropbox/MPC/Checker/kepcart/libkepcart.so")
lib = CDLL("/home/mpayne/Dropbox/MPC/MPChecker/kepcart/libkepcart.so")
   
def keplerian(GM, state):
    """
    Computes the Keplerian orbital elements a, e, incl, longnode,
    argperi, and meananom, given a GM constant and an input state.
    
    *Returns*
        (a, e, incl, longnode, argperi, meananom) : tuple of floats
    
    """

    _keplerian = lib.keplerian
    _keplerian.argtypes = (c_double, State)
    _keplerian.restype = None

    a = c_double()
    e = c_double()
    incl = c_double()
    longnode = c_double()
    argperi = c_double()
    meananom = c_double()

    return_value = _keplerian(GM, state, byref(a), byref(e), byref(incl), byref(longnode), byref(argperi), byref(meananom))

    return (a.value, e.value, incl.value, longnode.value, argperi.value, meananom.value)

def keplerians(num, GM, state_arr):
    """
    Computes arrays of Keplerian orbital elements a, e, incl, longnode,
    argperi, and meananom, given a GM constant and an array of input states.
    
    *Returns*
    numpy arrays of a, e, incl, longnode, argperi, meananom
    
    """

    StateArray = State * num

    a_arr = np.zeros((num), dtype=np.double)
    e_arr = np.zeros((num), dtype=np.double)
    incl_arr = np.zeros((num), dtype=np.double)
    longnode_arr = np.zeros((num), dtype=np.double)
    argperi_arr = np.zeros((num), dtype=np.double)
    meananom_arr =np.zeros((num), dtype=np.double)

    _keplerians = lib.keplerians
    _keplerians.argtypes = (c_int, c_double, POINTER(StateArray), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double))
    _keplerians.restype = None


    return_value = _keplerians(num, GM, byref(state_arr),
                               a_arr.ctypes.data_as(POINTER(c_double)),
                               e_arr.ctypes.data_as(POINTER(c_double)),
                               incl_arr.ctypes.data_as(POINTER(c_double)),
                               longnode_arr.ctypes.data_as(POINTER(c_double)),
                               argperi_arr.ctypes.data_as(POINTER(c_double)),
                               meananom_arr.ctypes.data_as(POINTER(c_double)))

    return a_arr, e_arr, incl_arr, longnode_arr, argperi_arr, meananom_arr




def cartesian(GM, a, e, incl, longnode, argperi, meananom):
    """
    Computes the cartesian state given a GM constant and the orbital elemments 
    a, e, incl, longnode, argperi, and meananom.
    
    *Returns*
        (x, y, z, xd, yd, zd) : tuple of floats
    
    """

    #lib = CDLL('/Users/mholman/slib_test/libkepcart.so')
    
    _cartesian = lib.cartesian
    _cartesian.argtypes = (c_double, c_double, c_double, c_double, c_double, c_double, c_double, POINTER(State))
    _cartesian.restype = None

    state = State(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    return_value = _cartesian(GM, a, e, incl, longnode, argperi, meananom, byref(state))

    return state

def cartesians(num, GM, a_arr, e_arr, incl_arr, longnode_arr, argperi_arr, meananom_arr):
    """
    Computes the cartesian states given a GM constant and arrays of the orbital elemments 
    a, e, incl, longnode, argperi, and meananom.
    
    *Returns*
    numpy arrays of x, y, z, xd, yd, zd states. : tuple of floats
    
    """

    StateArray = State * num
    state_arr = StateArray()

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


def cartesian_vectors(num, GM, a_arr, e_arr, incl_arr, longnode_arr, argperi_arr, meananom_arr):
    """
    Computes the cartesian position and velocity vectors given a GM constant and arrays of the orbital elemments 
    a, e, incl, longnode, argperi, and meananom.
    
    *Returns*
    arrays of x, y, z and xd, yd, zd values.
    
    """

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
                                      byref(vel_arr))

    return pos_arr, vel_arr


def cartesian_elements(num, GM, elements_arr):
    """
    Computes the cartesian position and velocity vectors given a GM constant and an array of sets of orbital elemments 
    a, e, incl, longnode, argperi, and meananom.
    
    *Returns*
    arrays of x, y, z and xd, yd, zd values.
    
    """

    ElementsArray = Elements * num
    
    size = num*3
    array_of_size_doubles = c_double*size

    pos_arr = array_of_size_doubles()
    vel_arr = array_of_size_doubles()

    _cartesian_elements = lib.cartesian_elements
    _cartesian_elements.argtypes = (c_int, c_double, POINTER(ElementsArray), POINTER(array_of_size_doubles), POINTER(array_of_size_doubles))
    _cartesian_elements.restype = None

    return_value = _cartesian_elements(num, GM, byref(elements_arr),
                                      byref(pos_arr),
                                      byref(vel_arr))

    return pos_arr, vel_arr


def universal_kepler(GM, state0, dt):
    """
    Takes state0 and advanced by time dt to state1.
    Assume a GM constant (mass of sun, for example)
    
    *Returns* 
        state1: State struct
    """

    _universal_kepler = lib.universal_kepler
    _universal_kepler.argtypes = (c_double, c_double, POINTER(State), POINTER(State))
    _universal_kepler.restype = None

    state1 = State(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    return_value = _universal_kepler(GM, dt, byref(state0), byref(state1))

    return state1
