# /mpcutilities/tests/test_classes_kepcart.py

"""

--------------------------------------------------------------

Oct 2018

Payne


Test the classes that are used by the kepcart conversion funcs


Tests are organized as follows

(i) Basic tests of class creation and method availability

(ii) ...

 Note that tests of kepcart conversions are in:
 /mpcutilities/tests/test_kepcart.py

 --------------------------------------------------------------

    """


# Import third-party packages
# --------------------------------------------------------------
#import unittest
import numpy as np

# Importing of local modules/packages required for this test
# --------------------------------------------------------------
import mpcutilities.phys_const as PHYS

# Import the specific package/module/function we are testing 
# --------------------------------------------------------------
import mpcutilities.classes as Classes


# (i) Basic test(s) of the classes required for the kepcart package
# ----------------------------------------------------------

def test_creation_of_cartstate():
    '''
    Create a single cartesian state
    '''
    # Define example cartesian coordinates
    x,y,z,xd,yd,zd = (4.545057765534996E-01  , 2.658978295316225E+00 , 1.275899778527264E-04, -1.040167152290219E-02  , 1.015042372911457E-03 , 1.948608211319801E-03)
    # Test that a state is being generated
    S      = Classes.CartState(x,y,z,xd,yd,zd)
    assert(  isinstance(S, Classes.CartState ) )
    assert( np.all(  np.array((S.x,S.y,S.z,S.xd,S.yd,S.zd))==np.array((x,y,z,xd,yd,zd))))
    
    
def test_standard_methods_cartstate():
    '''
    Test attributes / methods of CartState
    '''
    # Create a single cartesian state (as per "test_creation_of_cartstate()" above)
    x,y,z,xd,yd,zd = (4.545057765534996E-01  , 2.658978295316225E+00 , 1.275899778527264E-04, -1.040167152290219E-02  , 1.015042372911457E-03 , 1.948608211319801E-03)
    S      = Classes.CartState(x,y,z,xd,yd,zd)

    # Posn (xyz) method
    assert( np.all( np.array( [S.get_xyz()] ) == np.array( [x,y,z] ) ) )

    # Vel (uvw) method
    assert( np.all( np.array( [S.get_uvw()] ) == np.array( [xd,yd,zd] ) ) )

    # Rotation method
    R      = Classes.CartState(x,y,z,xd,yd,zd)
    R.ec_to_eq()
    assert( S.x == R.x ) 
    assert( S.y != R.y ) 
    assert( S.z != R.z ) 
    assert( np.all(  np.array((R.x,R.y,R.z)) ==  np.array(np.dot(PHYS.rot_mat_ec_to_eq , [S.x,S.y,S.z] )) ) )






def test_creation_of_cartstatearray():
    '''
    Test creation of an array of cartesian states
    '''
    # Create a single cartesian state (as per "test_creation_of_cartstate()" above)
    x,y,z,xd,yd,zd = (4.545057765534996E-01  , 2.658978295316225E+00 , 1.275899778527264E-04, -1.040167152290219E-02  , 1.015042372911457E-03 , 1.948608211319801E-03)
    S              = Classes.CartState(x,y,z,xd,yd,zd)

    # Now make an array of cartesian states of length 2
    N              = 2 
    StateArray     =  Classes.CartState * N
    SA             = (Classes.CartState * N)() # CartStateArray()
    SA[0]          = S

    assert(  isinstance(SA[0], Classes.CartState ) )
    assert( len(SA) == N)
    assert( np.all(  np.array((SA[0].x,SA[0].y,SA[0].z,SA[0].xd,SA[0].yd,SA[0].zd))==np.array((x,y,z,xd,yd,zd))))

    


def test_creation_of_kepstate():
    '''
    Test creation of a single KepState object
    '''
    # Define example keplerian elements
    a,e,i,O,o,M = (2.781870259390299, 0.07692192514945771, 0.18480538653567896, 1.4012469961929193, 1.237855673926063, -1.0950384781408455)
    KS          = Classes.KepState(a,e,i,O,o,M)
    # Test that a state is being generated
    assert(  isinstance(KS, Classes.KepState ) )
    assert(  np.all(  np.array((KS.a , KS.e , KS.incl , KS.longnode , KS.argperi , KS.meananom , ))==np.array((a,e,i,O,o,M))))
    

def test_creation_of_kepstatearray():
    '''
    Test creation of an array of keplerian states
    '''
    # Create an array of keplerian states of length 3
    a,e,i,O,o,M = (2.781870259390299, 0.07692192514945771, 0.18480538653567896, 1.4012469961929193, 1.237855673926063, -1.0950384781408455)
    N           = 3 
    KS          = (Classes.KepState * N)()
    for n in range(N):
        KS[n]   = Classes.KepState(a,e,i,O,o,M)
    
    # Test that the kepStateArray has the expected properties
    assert ( len(KS) ==  N)
    assert(  isinstance( KS[0], Classes.KepState ) )
    for n in range(N):
        assert(  np.all(
                            np.array((KS[n].a ,KS[n].e ,KS[n].incl ,KS[n].longnode ,KS[n].argperi ,KS[n].meananom , ))==
                            np.array((a,e,i,O,o,M))
                        )
               )



def test_creation_of_stateepoch():
    '''
    Test the creation of a combined CartStateEpoch structure to carry both the cartesian state & the epoch of validity

    '''
    # Create a single CartStateEpoch
    x,y,z,xd,yd,zd,JD = (4.545057765534996E-01, 2.439514137005412E+00, 1.057797885511556E+00, -1.040167152290219E-02, 1.561713370620346E-04, 2.191573748133731E-03, 2457933.5)
    CSE    = Classes.CartStateEpoch( Classes.CartState(x,y,z,xd,yd,zd), JD )
    # Check that the created CartStateEpoch, CSE, has the expected properties
    assert(  isinstance(CSE,       Classes.CartStateEpoch ) )
    assert(  isinstance(CSE.CartState, Classes.CartState ) )



def test_creation_of_stateepocharray():
    '''
    Test the creation of an array of combined CartStateEpoch structures
        
    '''
    # Create a single CartStateEpoch as in "test_creation_of_stateepoch" above
    x,y,z,xd,yd,zd,JD = (4.545057765534996E-01, 2.439514137005412E+00, 1.057797885511556E+00, -1.040167152290219E-02, 1.561713370620346E-04, 2.191573748133731E-03, 2457933.5)
    CSE    = Classes.CartStateEpoch( Classes.CartState(x,y,z,xd,yd,zd), JD )
    
    # Create a CartStateEpochArray of length 3 and populate it using the single CartStateEpoch, CSE
    N = 3
    SEA   = (Classes.CartStateEpoch * N)()
    for n in range(N): 
        SEA[n]=CSE
    
    # Checking StateEpochArray has the expected properties
    assert( len(SEA) == N)
    assert(  isinstance(SEA[0], Classes.CartStateEpoch ) )
    assert(  np.allclose( [SE.epoch for SE in SEA] , [ JD for n in range(N)] )  )    ### <<<--- WHY WOULD THIS WORK: JD != 2457933.5 ???
    for n in range(N):
        assert(
               np.all(
                        np.array((SEA[n].CartState.x , SEA[n].CartState.y , SEA[n].CartState.z , SEA[n].CartState.xd , SEA[n].CartState.yd , SEA[n].CartState.zd , SEA[n].epoch )) ==
                        np.array((x,y,z,xd,yd,zd,JD))
                      )
                )




def test_creation_of_KepStateEpoch():
    '''
    The creation of a combined KepStateEpoch structure to carry both the keplerian elements & the epoch of validity
    '''
    # Create single KepStateEpoch
    a,e,i,O,o,M,JD = (2.781870259390299, 0.07692192514945771, 0.18480538653567896, 1.4012469961929193, 1.237855673926063, -1.0950384781408455, 2457933.5)
    KSE    = Classes.KepStateEpoch( Classes.KepState(a,e,i,O,o,M), JD )
    # Checking single KepStateEpoch, KSE, has the expected properties
    assert(  isinstance(KSE,          Classes.KepStateEpoch ) )
    assert(  isinstance(KSE.KepState, Classes.KepState ) )




def test_creation_of_KepStateEpochArray():
    '''
    Test the creation of an array of combined KepStateEpoch structures
        
    '''
    # Create a single KepStateEpoch as in "test_creation_of_KepStateEpoch" above
    a,e,i,O,o,M,JD = (2.781870259390299, 0.07692192514945771, 0.18480538653567896, 1.4012469961929193, 1.237855673926063, -1.0950384781408455, 2457933.5)
    KSE    = Classes.KepStateEpoch( Classes.KepState(a,e,i,O,o,M), JD )

    # Create a KepStateEpochArray of length 3 and populate it using the single KepStateEpoch, KSE
    N = 3
    KSEA   = (Classes.KepStateEpoch * N)()
    for n in range(N):
        KSEA[n]=KSE
    
    # Checking KepStateEpochArray
    assert( len(KSEA) == N)
    assert(  isinstance(KSEA[0], Classes.KepStateEpoch ) )
    assert(  np.allclose( [KSE.epoch for KSE in KSEA] , [ JD for n in range(N)] )  )
    for n in range(N):
        assert(
               np.all(
                      np.array((KSEA[n].KepState.a , KSEA[n].KepState.e , KSEA[n].KepState.incl , KSEA[n].KepState.longnode , KSEA[n].KepState.argperi , KSEA[n].KepState.meananom , KSEA[n].epoch )) ==
                      np.array((a,e,i,O,o,M,JD))
                      )
               )






