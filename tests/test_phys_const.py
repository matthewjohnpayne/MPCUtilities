# /mpcutilities/tests/test_kepcart.py

"""
# --------------------------------------------------------------
# Oct 2018
# Payne
#
# Test the random assortedment of physical constants
# and conversion factors that I have collected in
# /mpcutilities/mpcutilities/phys_const.py
#
# --------------------------------------------------------------

"""

# Import third-party packages
# --------------------------------------------------------------
import numpy as np

# Import the specific package/module/function we are testing 
# --------------------------------------------------------------
import mpcutilities.phys_const  as PHYS

# Basic test(s) of the "phys_const" package
# ----------------------------------------------------------
def test_existance_of_constants():
    '''
    Test that the expected constants exist 
    '''
    assert(  isinstance(PHYS.GMsun,            float) )
    assert(  isinstance(PHYS.au_km,            float) )
    assert(  isinstance(PHYS.speed_of_light,   float) )
    assert(  isinstance(PHYS.ecl,              float) )
    assert(  isinstance(PHYS.Rearth_km,        float) )
    assert(  isinstance(PHYS.Rearth_AU,        float) )
    assert(  isinstance(PHYS.hr2deg,           float) )

    assert(  isinstance(PHYS.LLT_Tol,          float) )


def test_rotn_matricees_and_functions():
    '''
    Test that the rotations and matricees work as expected
    '''

    # Check Rotn Matricees exist
    assert(  isinstance( PHYS.rot_mat_ec_to_eq, np.ndarray ) )
    assert(  isinstance( PHYS.rot_mat_eq_to_ec, np.ndarray ) )
    assert(  isinstance( PHYS.rotate_matrix(PHYS.ecl), np.ndarray ) )
    
    # Check that the rotn matricees have the expected components ...
    R = [[ 1.,          0.,          0.        ],
         [ 0.,          0.91748213, -0.39777699],
         [ 0.,          0.39777699,  0.91748213]
         ]
    assert(  np.allclose( R , PHYS.rot_mat_ec_to_eq,  rtol=1e-09 ))

    # Check that the rotation matricees work "as expected"
    # I.e. operate on a vector and give a reasonable answer
    Cxyz = np.array([1,0,0])
    ExpectedQxyz = np.array([1,0,0])
    Qxyz = np.dot(PHYS.rot_mat_ec_to_eq , Cxyz )
    assert( np.allclose( ExpectedQxyz, Qxyz ) )

    # Check that the matricees give the same answer as the functions
    Cxyz = np.array([1,0,0])
    MQxyz = np.dot(PHYS.rot_mat_ec_to_eq , Cxyz )
    FQxyz = PHYS.rotate_ec_to_eq(Cxyz)
    assert( np.allclose( MQxyz, FQxyz ) )
    Qxyz = np.array([1,0,0])
    MCxyz = np.dot(PHYS.rot_mat_ec_to_eq , Qxyz )
    FCxyz = PHYS.rotate_ec_to_eq(Qxyz)
    assert( np.allclose( MCxyz, FCxyz ) )
    
    # Check that "rot_mat_ec_to_eq" & "rot_mat_eq_to_ec" cancel one another out
    Cxyz = np.array([1,0,0])
    Qxyz = np.dot(PHYS.rot_mat_ec_to_eq , Cxyz )
    outCxyz = np.dot(PHYS.rot_mat_eq_to_ec , Qxyz )
    assert( np.allclose( Cxyz, outCxyz ) )


