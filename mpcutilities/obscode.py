# /mpcutilities/mpcutilities/obscode.py
"""
    
--------------------------------------------------------------

Nov 2018

Payne

 - Original by Sonia Keys 

Notes by Payne:

 - This seems to be yet another set of functionalities by Sonia that are almost
 identical to functions written by Holman in mpc_library/mpcutilities
 
 - Needs to be rationalized

Original notes by Keys:

 - Function read5c loads file obscode.dat.  Currently supported is the current
style 5 column obscode.dat format, the format currently publicly accessible at
http://www.minorplanetcenter.net/iau/lists/ObsCodes.html.

 - The higher lever function siteXYZ is designed to take the output of read5c,
but should be general enough to also take output of a reader of the new
obscode.dat format.


--------------------------------------------------------------

"""

# Import third-party packages
# --------------------------------------------------------------
import collections
import math

__all__ = ['Site5c', 'read5c', 'siteXYZ']


Site5c = collections.namedtuple('Site5c', [
    'blank',  # bool:   True if the parallax data is completely blank
    'long',   # float:  parallax constant: east longitude of observatory in deg
    'rcos',   # float:  parallax constant: ρ cos φ′
    'rsin',   # float:  parallax constant: ρ sin φ′
    'site',   # string: site name
])
"""named tuple with fields for the five columns of obscode.dat

Parameters
----------
blank : bool
    True if the parallax data is completely blank
long : float
    parallax constant: east longitude of observatory in degrees
rcos : float
    parallax constant: ρ cos φ′
rsin : float
    parallax constant: ρ sin φ′
site : string
    site name
"""


def parseLine(line):
    """
    Parses a single line of obscode.dat.

    Parameters
    ----------
    line : str
        A line of at least 30 characters.  A line is valid if all three
        parallax constants (long, rcos, and rsin) parse or if all three are
        blank.

    Returns
    -------
    Site5c
        Site5c object for valid lines, None otherwise.
    """
    pc = line[4:30]
    blank = pc.isspace()
    if blank:
        long, rcos, rsin = 0., 0., 0.
    else:
        try:
            long = float(pc[:9])
            rcos = float(pc[9:17])
            rsin = float(pc[17:])
        except ValueError:
            return None  # invalid line

    return Site5c(
        blank=blank,
        long=long,
        rcos=rcos,
        rsin=rsin,
        site=line[30:].strip()
    )


def read5c(fn='obscode.dat'):
    """
    reads a file current style (5-column) obscode.dat into a dictionary.

    This is the format currently publicly accessible at
    http://www.minorplanetcenter.net/iau/lists/ObsCodes.html.

    Parameters
    ----------
    fn : string
        filename to read, default is obscode.dat

    Returns
    -------
    dictionary
        Keys are the 3-character obscode strings, values are Site5c named tuples
        of the parsed data.

    Examples
    --------
    >>> sites = read5c('obscode.dat')
    >>> sites['G96'].site
    'Mt. Lemmon Survey'
    """
    with open(fn) as f:
        sites = {}
        for line in f:
            s = parseLine(line)
            if s:
                sites[line[:3]] = s
        return sites


def siteXYZ(sites):
    """
    computes geocentric site coordinates from parallax constants.

    Parameters
    ----------
    sites : dictionary
        Typically the value returned from read5c.  Keys are obscodes and values
        must have attributes blank, rcos, rsin, and long.

    Returns
    -------
    dictionary
        Keys are the keys of the sites parameter, values are 3-tuples floats,
        the geocentric coordinates, where blank = False and the 3-tuple
        (None, None, None) where blank = True.
    """
    ObservatoryXYZ = {}
    for code, s in sites.items():
        if not s.blank:
            longitude = s.long * math.pi / 180
            x = s.rcos * math.cos(longitude)
            y = s.rcos * math.sin(longitude)
            z = s.rsin
            ObservatoryXYZ[code] = (x, y, z)
        else:
            ObservatoryXYZ[code] = (None, None, None)
    return ObservatoryXYZ


if __name__ == '__main__':
    import doctest
    doctest.testmod()
