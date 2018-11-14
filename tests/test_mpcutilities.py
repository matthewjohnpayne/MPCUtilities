# /mpcutilities/tests/test_mpcutilities.py

"""
    
--------------------------------------------------------------

Oct 2018

Payne

Test the classes & functions in the mpcutilities module

Lots of these functions are from Holman 
(and used to be in MPC_library)

*** MORE TESTS NEED TO BE CREATED ***

--------------------------------------------------------------

"""


# Import third-party packages
# --------------------------------------------------------------
#import unittest
import numpy as np
from pkg_resources import resource_exists, resource_filename

# Importing of local modules/packages required for this test
# --------------------------------------------------------------
import mpcutilities.phys_const as PHYS

# Import the specific package/module/function we are testing
# --------------------------------------------------------------
import mpcutilities.mpcutilities as MPCU


# (i) Basic test(s) of Observatory-Class
# ----------------------------------------------------------

def test_ObsCodesFile_exists_and_is_populated():
    '''
        ...
    '''
    # Check file exists
    assert resource_exists('mpcutilities','data/ObsCodes.txt')

    # Read the file
    with open(resource_filename('mpcutilities','data/ObsCodes.txt'), 'r') as f:
        data = f.readlines()
    
    # Check the  file is long enough
    assert len(data) > 2000

    # Check the data looks reasonable
    assert "Code  Long.   cos      sin    Name" in data[0]
    assert "000   0.0000 0.62411 +0.77873 Greenwich" in data[1]


def test_parseObsCode():
    '''
    ...
    '''
    O = MPCU.Observatory()
    testLine = "000   0.0000 0.62411 +0.77873 Greenwich"
    output = O.parseObsCode(testLine)
    assert len(output) == 5
    assert testLine[:3]    == output[0]
    assert testLine[4:13]  == output[1]
    assert testLine[13:21] == output[2]
    assert testLine[21:30] == output[3]
    assert testLine[30:].rstrip('\n') == output[4]


def test_ObservatoryClass_creation():
    '''
    Test the instantiation of an Observatory-class object
    
    On initialization, this opens a file and gets the geocentric-cartesian posns of each obscode
    '''
    # Create class (causes file-open & cache)
    # Test observatory dictionary got created
    O = MPCU.Observatory()
    assert isinstance( O, MPCU.Observatory )
    assert isinstance( O.ObservatoryXYZ , dict )

    # Test observatory dictionary has the expected obs-codes
    for obsCode in ['G96','T05','F51']:
        assert obsCode in O.ObservatoryXYZ.keys()
    assert len(O.ObservatoryXYZ) > 2000

    # Test observatory dictionary has the expected geocentric x,y,z positions for select obs-codes
    assert O.ObservatoryXYZ['000'] == (0.62411, 0.0, 0.77873)  ## <<-- Greenwich, => same as inputs


def test_parseXYZ():
    pass

def test_getObservatoryPosition():
    pass

def test_getSatellitePosition():
    pass


# (ii) Basic test(s) of EarthAndTime-Class
# ----------------------------------------------------------

def test_EarthAndTime_creation():
    pass

def test_polarmotion_functions():
    pass

def test_timeTransformations():
    pass

def test_getEarthPosition():
    pass

def test_dateFormatConversions():
    pass

def test_MagnitudeConversions():
    pass

