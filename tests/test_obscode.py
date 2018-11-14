# mpcutilities/tests/test_obscode.py
"""
    
--------------------------------------------------------------

Nov 2018

Payne

- Original by Sonia Keys

Notes by Payne:

- This seems to be yet another set of functionalities by Sonia that are almost
identical to functions written by Holman in mpc_library/mpcutilities

- Needs to be rationalized

--------------------------------------------------------------

"""




# Import third-party packages
# --------------------------------------------------------------
import novas.compat as novas
import pytest
from pkg_resources import resource_filename
import numpy as np

# Import neighboring packages
# --------------------------------------------------------------
import mpcutilities.leapsec as leapsec

# Import the specific package/module/function we are testing
# --------------------------------------------------------------
import mpcutilities.obscode as obscode

def test_leap():
    """example from leap-seconds.list"""
    l = leapsec.LeapSeconds()
    jd9 = novas.julian_date(1972, 6, 30, 23 + (3599. / 3600))
    assert l.getLeapSeconds(jd9) == 10
    jd0 = novas.julian_date(1972, 7, 1, 0)
    assert l.getLeapSeconds(jd0) == 11


@pytest.fixture(scope='module')
def sites():
    # look for a real obcode file
    ObsFilename = resource_filename('mpcutilities','data/ObsCodes.txt')
    for fn in [ObsFilename,
            '/var/www/html/mpc/public/iau/lists/ObsCodes.html']:
        try:
            return obscode.read5c(fn)
        except IOError:
            pass
    sys.exit("Could not open ObsCode file ")


def test_sites(sites):
    assert sites['G96'] == obscode.Site5c(blank=False, long=249.21128,
        rcos=0.845111, rsin=0.533614, site='Mt. Lemmon Survey')
    # Blank parallax constants
    assert sites['248'] == obscode.Site5c(blank=False, long=0.0, rcos=0.0,
        rsin=0.0, site='Hipparcos')
    # 0 parallax constants (not blank)
    assert sites['248'] == obscode.Site5c(blank=False, long=0.0, rcos=0.0,
        rsin=0.0, site='Hipparcos')
    # Southern hemisphere
    assert sites['413'] == obscode.Site5c(blank=False, long=149.06608,
        rcos=0.855595, rsin=-0.516262, site='Siding Spring Observatory')


def test_xyz(sites):
    xyz = obscode.siteXYZ(sites)
    assert np.allclose(xyz['W14'] , (0.04652656088240064, -0.8205960653868969, 0.567774) )
    assert np.allclose(xyz['248'] , (0.0, 0.0, 0.0) )
    assert np.allclose(xyz['413'] , (-0.7338957925529292, 0.43981788242192676, -0.516262) )
    assert xyz['250'] == (None, None, None)
