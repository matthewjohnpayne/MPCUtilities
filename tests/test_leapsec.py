# mpcutilities/tests/test_leapsec.py
"""
    
--------------------------------------------------------------

Nov 2018

Payne

Original code by Sonia Keys

Tests leapsec.py code by Sonia Keys
 - moved to /mpcutilities/mpcutilities/UNUSED_leapsec.py

*** Payne: I believe this code to be unnecessary / unused ***

*** Payne: I am keeping it for archival purposes ***

--------------------------------------------------------------

"""

# Third-party imports
import novas.compat as novas

# Import the specific package/module/function we are testing
import mpcutilities.leapsec as leapsec

def test_leap():
    """example from leap-seconds.list"""
    l = leapsec.LeapSeconds()
    jd9 = novas.julian_date(1972, 6, 30, 23 + (3599. / 3600))
    assert l.getLeapSeconds(jd9) == 10
    jd0 = novas.julian_date(1972, 7, 1, 0)
    assert l.getLeapSeconds(jd0) == 11
