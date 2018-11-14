# /mpcutilities/mpcutilities/leapsec.py
"""
    
--------------------------------------------------------------

Nov 2018

Payne

Original code by Sonia Keys 

Provides leap second information by year.  It requires a file
of leap second data in the format prepared by Judah Levine at NIST.
See https://www.ietf.org/timezones/data/leap-seconds.list for example.

*** Payne: I believe this code to be unnecessary / unused ***

*** Payne: I am keeping it for archival purposes ***

--------------------------------------------------------------

"""

# Import third-party packages
# --------------------------------------------------------------
import numpy as np
from pkg_resources import resource_filename

__all__ = ["LeapSeconds"]


class LeapSeconds:
    """loads leap second data from a file.

    Parameters
    ----------
    fn : string
        File name of leap second data, by default leap-seconds.list.
    jd_ref_utc : float
        Reference date, JD at which to start counting leap seconds.
        The default value of 2415020.5 is year 1900.0.
    """

    def __init__(self, fn=resource_filename('mpcutilities','data/leap-seconds.list'), jd_ref_utc=2415020.5):
        # Create two numpy arrays, secSinceArray and leapSecArray,
        # that I can use to make a look-up table.
        # This is just reading Judah Levine's file.
        self.jd_ref_utc = jd_ref_utc
        secSinceList = []
        leapSecList = []
        with open(fn) as f:
            for line in f:
                if not line.startswith('#'):
                    secondsSince, leapSeconds = line.rsplit('#')[0].split()
                    secSinceList.append(int(secondsSince))
                    leapSecList.append(int(leapSeconds))
        self.secSinceArray = np.array(secSinceList)
        self.leapSecArray = np.array(leapSecList)

    def getLeapSeconds(self, jd_utc):
        """
        gets number of leap seconds since the reference Julian date.

        Parameters
        ----------
        jd_utc : float
            JD of end of period since reference date.

        Returns
        -------
        int
            number of leap seconds since reference date.
        """
        # Given a Julian Date in UTC, determine the number of seconds
        # that have elapsed since the reference time.  Then find the
        # correspoding number of leap seconds.
        max_idx = self.secSinceArray.shape[0] - 1
        secsSince = self.secondsSinceRef(jd_utc)
        idx = np.searchsorted(self.secSinceArray, int(secsSince),
            side='right') - 1
        if idx < 0:
            return 0
        if idx > max_idx:
            return self.leapSecArray[max_idx]
        else:
            return self.leapSecArray[idx]

    def secondsSinceRef(self, jd_utc):
        """Determines the seconds (not just leap seconds) that have elapsed
        since the reference Julian date.

        Parameters
        ----------
        jd_utc : float
            JD of end of period since reference date.

        Returns
        -------
        float
            number of seconds since reference date.
        """
        secondsSince = (jd_utc - self.jd_ref_utc) * 24.0 * 60 * 60
        return secondsSince
