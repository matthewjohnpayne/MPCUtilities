# /mpcutilities/tests/test_kepcart.py

"""
# --------------------------------------------------------------
# Oct 2018
# Payne
#
# Test the kepcart conversion funcs that are in
# /mpcutilities/mpcutilities/kepcart.py
# 
# Note that kepcart needs classes that are defined in 
# /mpcutilities/mpcutilities/classes.py
# The classes are tested separately in:
# /mpcutilities/tests/test_classes_kepcart.py
#
# *** MORE TESTING OF THE ACCURACY OF THE CODE IS REQUIRED ***
#    ************* ESP. FOR HYPERBOLIC CASES **************
# -- there is a copy of xv2els in /Users/matthewjohnpayne/Dropbox/kepcart.v2/ that may be of use for comparison-tests
#
# --------------------------------------------------------------

"""


# Import third-party packages
# --------------------------------------------------------------
#import unittest
import numpy as np

# Importing of local modules/packages required for this test
# --------------------------------------------------------------
import mpcutilities.phys_const  as PHYS
import mpcutilities.classes     as Classes


# Import the specific package/module/function we are testing
# --------------------------------------------------------------
import mpcutilities.kepcart as kc


# Basic test(s) of the kepcart package
# ----------------------------------------------------------

def test_cart2kep():
    '''
    Test the accuracy of the cart2kep function
    - This converts a CartState to individual Keplerian elements
    '''
    # Sample data for Ceres from JPL
    S      = Classes.CartState(4.54505776553499E-01, 2.65897829531622E+00, 1.27589977852726E-04, -1.04016715229021E-02, 1.01504237291145E-03, 1.94860821131980E-03)
    aeiOoM = (2.781870259390299, 0.07692192514945771, 0.18480538653567896, 1.4012469961929193, 1.237855673926063, -1.0950384781408455)

    # Use cart2kep to calculate the elements
    els    = kc.cart2kep(PHYS.GMsun , S)
    
    # Test datatype & accuracy of returned elements
    assert(  isinstance( els,  tuple ) )
    assert( len(els) == 6)
    assert( np.allclose( els , aeiOoM )) 

def test_keplerian():
    '''
    keplerian should be exactly the same as "cart2kep"
     - test that it is
    (Provided for backward compatibility with Holman's legacy code)
    '''
    # Sample data for Ceres from JPL
    S      = Classes.CartState(4.54505776553499E-01, 2.65897829531622E+00, 1.27589977852726E-04, -1.04016715229021E-02, 1.01504237291145E-03, 1.94860821131980E-03)
    # Use cart2kep to calculate the elements
    elsC2K    = kc.cart2kep(PHYS.GMsun , S)
    # Use keplerian to calculate the elements
    elsKEP    = kc.keplerian(PHYS.GMsun , S)
    # Test that the results are the same
    assert np.allclose( elsC2K , elsKEP )





def test_cart2kep_array():
    '''
    Test the accuracy of the cart2kep_array function:
     - This converts an array of CartStates to arrays of individual Keplerian elements
    '''
    # Create array of CartStates of length 2
    # Sample data for Ceres from JPL
    N                 = 2
    SA                = (Classes.CartState * N)()
    x,y,z,xd,yd,zd    = (4.54505776553499E-01, 2.65897829531622E+00, 1.27589977852726E-04, -1.04016715229021E-02, 1.01504237291145E-03, 1.94860821131980E-03)
    aeiOoM            = (2.781870259390299, 0.07692192514945771, 0.18480538653567896, 1.4012469961929193, 1.237855673926063, -1.0950384781408455)
    SA[0]             = Classes.CartState(x,y,z,xd,yd,zd)
    SA[1]             = Classes.CartState(2.*x,2.*y,2.*z, 0.5*xd, 0.5*yd, 0.5*zd)
    
    # Use cart2kep_array to calculate the elements
    arrays_ = kc.cart2kep_array(PHYS.GMsun , SA)
    
    # Test datatype & accuracy of returned elements
    assert( isinstance(arrays_, tuple) )
    assert(       len(arrays_) == 6)
    a_,e_,i_,O_,o_,M_ = arrays_
    assert( isinstance(a_, np.ndarray) )
    assert( len(a_) == N)
    assert( np.allclose( (a_[0],e_[0],i_[0],O_[0],o_[0],M_[0]) , aeiOoM )) 


def test_keplerians():
    '''
        keplerians should be exactly the same as "cart2kep_array"
        - test that it is
        (Provided for backward compatibility with Holman's legacy code)
        '''
    # Create array of CartStates of length 2
    # Sample data for Ceres from JPL
    N                 = 2
    SA                = (Classes.CartState * N)()
    x,y,z,xd,yd,zd    = (4.54505776553499E-01, 2.65897829531622E+00, 1.27589977852726E-04, -1.04016715229021E-02, 1.01504237291145E-03, 1.94860821131980E-03)
    SA[0]             = Classes.CartState(x,y,z,xd,yd,zd)
    SA[1]             = Classes.CartState(2.*x,2.*y,2.*z, 0.5*xd, 0.5*yd, 0.5*zd)
    # Use cart2kep_array to calculate the elements
    arrays_C2K = kc.cart2kep_array(PHYS.GMsun , SA)
    # Use keplerians to calculate the elements
    arrays_KEP = kc.keplerians(PHYS.GMsun , SA)
    # Test that the results are the same
    assert np.allclose( arrays_C2K, arrays_KEP)
    # Just double-checking the test-comparison ...
    a_C2K,e_C2K,i_C2K,O_C2K,o_C2K,M_C2K = arrays_C2K
    a_KEP,e_KEP,i_KEP,O_KEP,o_KEP,M_KEP = arrays_KEP
    assert( np.allclose( (a_C2K,e_C2K,i_C2K,O_C2K,o_C2K,M_C2K) , (a_KEP,e_KEP,i_KEP,O_KEP,o_KEP,M_KEP) ))






def test_kep2cartState():
    '''
    Test the conversion of keplerian elements to cartesian state using "kep2cartState"
    '''
    
    # Create elements (sample data for Ceres from JPL Horizons)
    a,e,i,O,o,M       = (2.781870259390299, 0.07692192514945771, 0.18480538653567896, 1.4012469961929193, 1.237855673926063, -1.0950384781408455)
    x,y,z,xd,yd,zd    = (4.54505776553499E-01, 2.65897829531622E+00, 1.27589977852726E-04, -1.04016715229021E-02, 1.01504237291145E-03, 1.94860821131980E-03)
    
    # Use kep2cartState to calculate the cartesian elements
    state             = kc.kep2cartState(PHYS.GMsun , a,e,i,O,o,M)

    # Test datatype & accuracy of returned CartState
    assert(  isinstance( state,  Classes.CartState ) )
    assert( np.allclose(  np.array((state.x,state.y,state.z,state.xd,state.yd,state.zd)) , np.array((x,y,z,xd,yd,zd))))

def test_cartesian():
    '''
    cartesian should be exactly the same as "kep2cartState"
    - test that it is
    (Provided for backward compatibility with Holman's legacy code)
    '''
    # Create elements (sample data for Ceres from JPL Horizons)
    a,e,i,O,o,M = (2.781870259390299, 0.07692192514945771, 0.18480538653567896, 1.4012469961929193, 1.237855673926063, -1.0950384781408455)
    # Use kep2cartState to calculate the elements
    stateK2C    = kc.kep2cartState(PHYS.GMsun , a,e,i,O,o,M)
    # Use cartesian to calculate the elements
    stateCART   = kc.cartesian(PHYS.GMsun , a,e,i,O,o,M)
    # Test that the results are the same
    np.allclose(    np.array((stateK2C.x,stateK2C.y,stateK2C.z,stateK2C.xd,stateK2C.yd,stateK2C.zd)) ,
                    np.array((stateCART.x,stateCART.y,stateCART.z,stateCART.xd,stateCART.yd,stateCART.zd))
                )


def test_kep2cartStateArray():
    '''
    Test the conversion of an array of keplerian elements to cartesian state using "kep2cartStateArray"
    '''
    # Create elements (sample data for Ceres from JPL Horizons)
    x,y,z,xd,yd,zd    = (4.54505776553499E-01, 2.65897829531622E+00, 1.27589977852726E-04, -1.04016715229021E-02, 1.01504237291145E-03, 1.94860821131980E-03)
    a,e,i,O,o,M       = (2.781870259390299, 0.07692192514945771, 0.18480538653567896, 1.4012469961929193, 1.237855673926063, -1.0950384781408455)
    N  = 10
    a_ = np.full(N,a)
    e_ = np.full(N,e)
    i_ = np.full(N,i)
    O_ = np.full(N,O)
    o_ = np.full(N,o)
    M_ = np.full(N,M)
    
    # Use kep2cartStateArray to calculate the elements
    SA = kc.kep2cartStateArray(PHYS.GMsun , a_,e_,i_,O_,o_,M_)
    
    # Test datatype & accuracy of returned CartState
    assert(  isinstance( SA[0],  Classes.CartState ) )
    assert( len(SA) == N)
    for n in range(N):
        assert( np.allclose(  np.array((SA[n].x,SA[n].y,SA[n].z,SA[n].xd,SA[n].yd,SA[n].zd)) , np.array((x,y,z,xd,yd,zd))))

def test_cartesians():
    '''
    cartesians should be exactly the same as "kep2cartStateArray"
    - test that it is
    (Provided for backward compatibility with Holman's legacy code)
    '''
    # Create elements (sample data for Ceres from JPL Horizons)
    a,e,i,O,o,M       = (2.781870259390299, 0.07692192514945771, 0.18480538653567896, 1.4012469961929193, 1.237855673926063, -1.0950384781408455)
    N  = 10
    a_ = np.full(N,a)
    e_ = np.full(N,e)
    i_ = np.full(N,i)
    O_ = np.full(N,O)
    o_ = np.full(N,o)
    M_ = np.full(N,M)
    # Use kep2cartStateArray to calculate the elements
    SA_K2C  = kc.kep2cartStateArray(PHYS.GMsun , a_,e_,i_,O_,o_,M_)
    # Use cartesians to calculate the elements
    SA_CART = kc.cartesians(PHYS.GMsun , a_,e_,i_,O_,o_,M_)
    # Test that the results are the same
    for n in range(N):
        assert( np.allclose(
                            np.array((SA_K2C[n].x,SA_K2C[n].y,SA_K2C[n].z,SA_K2C[n].xd,SA_K2C[n].yd,SA_K2C[n].zd)) ,
                            np.array((SA_CART[n].x,SA_CART[n].y,SA_CART[n].z,SA_CART[n].xd,SA_CART[n].yd,SA_CART[n].zd))
                            ))






def test_kep2cartPV():
    '''
    Test the conversion of arrays of keplerian elements to cartesian position and velocity arrays using "kep2cartPV"
    '''
    # Create elements (sample data for Ceres from JPL Horizons)
    x,y,z,xd,yd,zd    = (4.54505776553499E-01, 2.65897829531622E+00, 1.27589977852726E-04, -1.04016715229021E-02, 1.01504237291145E-03, 1.94860821131980E-03)
    a,e,i,O,o,M       = (2.781870259390299, 0.07692192514945771, 0.18480538653567896, 1.4012469961929193, 1.237855673926063, -1.0950384781408455)
    N  = 10
    a_ = np.full(N,a)
    e_ = np.full(N,e)
    i_ = np.full(N,i)
    O_ = np.full(N,O)
    o_ = np.full(N,o)
    M_ = np.full(N,M)
    
    # Use kep2cartPV to calculate the posn & velocity arrays
    XYZ, UVW = kc.kep2cartPV(PHYS.GMsun , a_,e_,i_,O_,o_,M_)
    
    # Test datatype & accuracy of returned cartesian posn & velocity vectors
    assert(          len(XYZ) == N )
    assert(          len(UVW) == N )
    for xyz, uvw in zip ( XYZ, UVW ):
        assert( np.allclose( (xyz[0],xyz[1],xyz[2],uvw[0],uvw[1],uvw[2]) , (x,y,z,xd,yd,zd) ))

def test_cartesian_vectors():
    '''
    cartesian_vectors should be exactly the same as "kep2cartPV"
    - test that it is
    (Provided for backward compatibility with Holman's legacy code)
    '''
    # Create elements (sample data for Ceres from JPL Horizons)
    a,e,i,O,o,M       = (2.781870259390299, 0.07692192514945771, 0.18480538653567896, 1.4012469961929193, 1.237855673926063, -1.0950384781408455)
    N  = 10
    a_ = np.full(N,a)
    e_ = np.full(N,e)
    i_ = np.full(N,i)
    O_ = np.full(N,O)
    o_ = np.full(N,o)
    M_ = np.full(N,M)
    
    # Use kep2cartPV to calculate the posn & velocity arrays
    XYZ_K2C, UVW_K2C = kc.kep2cartPV(PHYS.GMsun , a_,e_,i_,O_,o_,M_)
    # Use cartesian_vectors to calculate the posn & velocity arrays
    XYZ_CART, UVW_CART = kc.cartesian_vectors(PHYS.GMsun , a_,e_,i_,O_,o_,M_)
    # Test that the results are the same
    for xyz_K2C, uvw_K2C,xyz_CART, uvw_CART  in zip ( XYZ_K2C, UVW_K2C, XYZ_CART, UVW_CART ):
        assert(
               np.allclose(
                        (xyz_K2C[0], xyz_K2C[1], xyz_K2C[2], uvw_K2C[0], uvw_K2C[1], uvw_K2C[2]),
                        (xyz_CART[0],xyz_CART[1],xyz_CART[2],uvw_CART[0],uvw_CART[1],uvw_CART[2])
                           )
               )




def test_kepState2cartPV():
    '''
    Test the conversion of an input array of KepState objects into cartesian position and velocity arrays using "kepState2cartPV"
    '''
    # Create elements (sample data for Ceres from JPL Horizons)
    a,e,i,O,o,M = (2.781870259390299, 0.07692192514945771, 0.18480538653567896, 1.4012469961929193, 1.237855673926063, -1.0950384781408455)
    N           = 3
    KSA          = (Classes.KepState * N)()
    for n in range(N):
        KSA[n]   = Classes.KepState(a,e,i,O,o,M)

    # Use kepState2cartPV to calculate the posn & velocity arrays
    XYZ, UVW = kc.kepState2cartPV( PHYS.GMsun , KSA)

    # Test datatype & accuracy of returned cartesian posn & velocity vectors
    assert( len(KSA) == N)
    assert(  isinstance(KSA[0], Classes.KepState ) )
    for n in range(N):
        # Checking the input is as expected
        assert(  np.all(
                        np.array((KSA[n].a , KSA[n].e , KSA[n].incl , KSA[n].longnode , KSA[n].argperi , KSA[n].meananom , ))==
                        np.array((a,e,i,O,o,M))
                        )
               )
        # Checking the output is as expected
        assert np.allclose(XYZ[n] , np.array([4.54505777e-01 ,2.65897830e+00, 1.27589978e-04]) )


def test_cartesian_elements():
    '''
    cartesian_elements should be exactly the same as "kepState2cartPV"
    - test that it is
    (Provided for backward compatibility with Holman's legacy code)
    # Create elements (sample data for Ceres from JPL Horizons)
    a,e,i,O,o,M = (2.781870259390299, 0.07692192514945771, 0.18480538653567896, 1.4012469961929193, 1.237855673926063, -1.0950384781408455)
    N           = 3
    KSA          = (Classes.KepState * N)()
    for n in range(N):
        KSA[n]   = Classes.KepState(a,e,i,O,o,M)
    
    # Use kepState2cartPV to calculate the posn & velocity arrays
    XYZ1, UVW1 = kc.kepState2cartPV( PHYS.GMsun , KSA)

    # Use cartesian_elements to calculate the posn & velocity arrays
    XYZ2, UVW2 = kc.kepState2cartPV( PHYS.GMsun , KSA)

    # Test that the returned values are the same
    for n in range(N):
        assert(  np.all(
                    np.array((KSA[n].a , KSA[n].e , KSA[n].incl , KSA[n].longnode , KSA[n].argperi , KSA[n].meananom , ))==
                    np.array((a,e,i,O,o,M))
                    )
               )
        '''
    pass

