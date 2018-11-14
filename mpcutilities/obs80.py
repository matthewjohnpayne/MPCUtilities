#!/usr/bin/env python3
# /mpcutilities/mpcutilities/obs80.py
"""

--------------------------------------------------------------

Oct 2018

Relocated by Payne

Original by Keys

Parses a file of 80 column observations line by line.


Returns list of named tuples from some short list of tuple types like
Optical, Satellite, and Radar.

Some string fields are parsed to numeric types, but otherwise
the tuple format corresponds directly to the 80 column format.

Significant overlap with Holman's MPC_library
 - Needs to be rationalized

Functions such as the following (written by Holman)
   might more naturally belong in here (obs80)

jdutc=date2JD(line1[15:32]),

ra=RA2hrRA(line1[32:44]),

dec=Dec2degDec(line1[44:56]),

--------------------------------------------------------------

"""

# Import third-party packages
# --------------------------------------------------------------
import collections
import math
import sys
import novas.compat as novas


# Import neighboring packages
# --------------------------------------------------------------
import mpcutilities.mpcutilities    as MPC
import mpcutilities.phys_const      as PHYS


# Methods for reading obs80-style observations
# --------------------------------------------------------------
__all__ = ['parse80', 'Optical', 'SpaceBased', 'Roving', 'Radar']


Optical = collections.namedtuple('Optical', [
    'num',    # string: Columns 1-5, packed perm desig + comet/sat stuff
    'desig',  # string: Columns 6-12, packed provisional or temporary desig
    'disc',   # string: Column 13, discovery asterisk
    'note1',  # string: Column 14, overloaded note/program code
    'note2',  # string: Column 15, overloaded record type/measurement tech
    'jdutc',  # float:  Columns 16-32, time of obs parsed to JD
    'ra',     # float:  Columns 33-44, RA parsed to hours
    'dec',    # float:  Columns 45-56, dec parsed to degrees
    'mag',    # float:  Columns 66-70, mag parsed to float
    'band',   # string: Column 71, mag band
    'cod'     # string: Columns 78-80, 3 character obs code.
])
"""named tuple type for the optical format documented at
http://www.minorplanetcenter.net/iau/info/OpticalObs.html.

Parameters
----------
num : string
    Columns 1-5, packed perm desig + comet/sat stuff
desig : string
    Columns 6-12, packed provisional or temporary desig
disc : string
    Column 13, discovery asterisk
note1 : string
    Column 14, overloaded note/program code
note2 : string
    Column 15, overloaded record type/measurement technology
jdutc : float
    Columns 16-32, time of obs parsed to JD
ra : float
    Columns 33-44, RA parsed to hours
dec : float
    Columns 45-56, dec parsed to degrees
mag : float
    Columns 66-70, mag parsed to float
band : string
    Column 71, mag band
cod : string
    Columns 78-80, 3 character obs code.
"""


SpaceBased = collections.namedtuple('SpaceBased', [
    # first 80 column record same as Optical except c15 is 'S'
    'num',    # string: Columns 1-5, packed perm desig + comet/sat stuff
    'desig',  # string: Columns 6-12, packed provisional or temporary desig
    'disc',   # string: Column 13, discovery asterisk
    'note1',  # string: Column 14, overloaded note/program code
    'jdutc',  # float:  Columns 16-32, time of obs parsed to JD
    'ra',     # float:  Columns 33-44, RA parsed to hours
    'dec',    # float:  Columns 45-56, dec parsed to degrees
    'mag',    # float:  Columns 66-70, mag parsed to float
    'band',   # string: Column 71, mag band
    'cod',    # string: Columns 78-80, 3 character obs code.
    # second 80 column record:
    'x',      # float:  Columns 35-45, converted to AU as needed
    'y',      # float:  Columns 47-57, converted to AU as needed
    'z'       # float:  Columns 59-69, converted to AU as needed
])
"""named tuple type for the space-based format documented at
http://www.minorplanetcenter.net/iau/info/SatelliteObs.html.

First 80 column record same as Optical except c15 note2 is omitted.
Second 80 column record contains spacecraft xyz coordinates.

Parameters
----------
num : string
    Columns 1-5, packed perm desig + comet/sat stuff
desig : string
    Columns 6-12, packed provisional or temporary desig
disc : string
    Column 13, discovery asterisk
note1 : string
    Column 14, overloaded note/program code
jdutc : float
    Columns 16-32, time of obs parsed to JD
ra : float
    Columns 33-44, RA parsed to hours
dec : float
    Columns 45-56, dec parsed to degrees
mag : float
    Columns 66-70, mag parsed to float
band : string
    Column 71, mag band
cod : string
    Columns 78-80, 3 character obs code.
x : float
    Columns 35-45 of second record, converted to AU as needed
y : float
    Columns 47-57 of second record, converted to AU as needed
z : float
    Columns 59-69 of second record, converted to AU as needed
"""


Roving = collections.namedtuple('Roving', [
    # first 80 column record same as Optical except c15 is 'V'
    'num',    # string: Columns 1-5, packed perm desig + comet/sat stuff
    'desig',  # string: Columns 6-12, packed provisional or temporary desig
    'disc',   # string: Column 13, discovery asterisk
    'note1',  # string: Column 14, overloaded note/program code
    'jdutc',  # float:  Columns 16-32, time of obs parsed to JD
    'ra',     # float:  Columns 33-44, RA parsed to hours
    'dec',    # float:  Columns 45-56, dec parsed to degrees
    'mag',    # float:  Columns 66-70, mag parsed to float
    'band',   # string: Column 71, mag band
    # second 80 column record:
    'lon',    # float:  Columns 35-44, longitude, parsed as degrees
    'lat',    # float:  Columns 46-55, latitude, parsed as degrees
    'alt'     # float:  Columns 57-61, altitude, parsed as meters
])
"""named tuple type for the roving observer format documented at
http://www.minorplanetcenter.net/iau/info/RovingObs.html.

First 80 column record same as Optical except c15 note2 is omitted.
Second 80 column record contains observer coordinates.

Parameters
----------
num : string
    Columns 1-5, packed perm desig + comet/sat stuff
desig : string
    Columns 6-12, packed provisional or temporary desig
disc : string
    Column 13, discovery asterisk
note1 : string
    Column 14, overloaded note/program code
jdutc : float
    Columns 16-32, time of obs parsed to JD
ra : float
    Columns 33-44, RA parsed to hours
dec : float
    Columns 45-56, dec parsed to degrees
mag : float
    Columns 66-70, mag parsed to float
band : string
    Column 71, mag band
lon : float
    Columns 35-44 of second record, longitude parsed as degrees
lat : float
    Columns 46-55 of second record, latitude parsed as degrees
alt : float
    Columns 57-61 of second record, altitude parsed as meters
"""


Radar = collections.namedtuple('Radar', [
    # first 80 column record same as Optical except c13 disc and c15 note2
    # are omitted.
    'num',      # string: Columns 1-5, packed perm desig + comet/sat stuff
    'desig',    # string: Columns 6-12, packed provisional or temporary desig
    'note1',    # string: Column 14, overloaded note/program code
    'jdutc',    # float:  Columns 16-32, time of obs parsed to JD
    'delay',    # float:  Columns 33-47, time delay in microsecs
    'doppler',  # float:  Columns 48-62, doppler shift in Hz
    'freq',     # float:  Columns 63-68 + 63-68 of second rec, tx frequency.
    'txcod',    # string: Columns 69-71, transmitter 3 character obs code.
    'rxcod',    # string: Columns 78-80, receiver 3 character obs code.
    # second 80 column record:
    'ret',       # string: Column 33, "S" for surface, "C" for mass center
    'udelay',    # float:  Columns 34-47, delay uncertainty in microsecs
    'udoppler',  # float:  Columns 48-62, doppler uncertainty in Hz
])
"""named tuple type for the radar observation format documented at
http://www.minorplanetcenter.net/iau/info/RadarObs.html.

Parameters
----------
num : string
    Columns 1-5, packed perm desig + comet/sat stuff
desig : string
    Columns 6-12, packed provisional or temporary desig
note1 : string
    Column 14, overloaded note/program code
jdutc : float
    Columns 16-32, time of obs parsed to JD
delay : float
    Columns 33-47, time delay in microsecs
doppler : float
    Columns 48-62, doppler shift in Hz
freq : float
    Columns 63-63 + 63-68 of second rec, tx frequency
txcod : string
    Columns 69-71, transmitter 3 character obs code
rxcod : string
    Columns 78-80, receiver 3 character obs code
ret : string
    Column 33 of second record, "S" for surface, "C" for mass center
udelay : float
    Columns 34-47 of second record, delay uncertainty in microsecs
udoppler : float
    Columns 48-62 of second record, doppler uncertainty in Hz
"""


def parseOpt(line):
    """ returns Optical for successful parse, raises exception otherwise"""
    n2 = note2(line)
    if n2 not in 'PeCTMcEOHNn AXx':
        raise ValueError('note2 = "{}"'.format(n2))

    return Optical(
        num=line[:5].strip(),
        desig=line[5:12].strip(),
        disc=line[12].strip(),
        note1=line[13].strip(),
        note2=n2.strip(),
        jdutc=date2JD(line[15:32]),
        ra=RA2hrRA(line[32:44]),
        dec=Dec2degDec(line[44:56]),
        mag=floatMag(line[65:70]),
        band=line[70].strip(),
        cod=line[77:80]
    )


def parseSat(line1, line2):
    """ returns Spacebased for successful parse, raises exception otherwise"""

    # if ... fields don't match
    #   raise ValueError("sat lines don't match")

    ax = floatPx(line2[34:45])
    ay = floatPx(line2[46:57])
    az = floatPx(line2[58:69])
    if line2[32] == '1':
        ax /= PHYS.au_km
        ay /= PHYS.au_km
        az /= PHYS.au_km

    return SpaceBased(
        # line1
        num=line1[:5].strip(),
        desig=line1[5:12].strip(),
        disc=line1[12].strip(),
        note1=line1[13].strip(),
        jdutc=date2JD(line1[15:32]),
        ra=RA2hrRA(line1[32:44]),
        dec=Dec2degDec(line1[44:56]),
        mag=floatMag(line1[65:70]),
        band=line1[70].strip(),
        cod=line1[77:80],
        # line2
        x=ax,
        y=ay,
        z=az
    )


def parseRov(line1, line2):
    """ returns Roving for successful parse, raises exception otherwise"""

    # if ... fields don't match
    #   raise ValueError("rov lines don't match")

    return Roving(
        # line1
        num=line1[:5].strip(),
        desig=line1[5:12].strip(),
        disc=line1[12].strip(),
        note1=line1[13].strip(),
        jdutc=date2JD(line1[15:32]),
        ra=RA2hrRA(line1[32:44]),
        dec=Dec2degDec(line1[44:56]),
        mag=floatMag(line1[65:70]),
        band=line1[70].strip(),
        # line2
        lon=float(line2[34:44]),
        lat=floatPx(line2[45:55]),
        alt=float(line2[56:61])
    )


def parseRad(line1, line2):
    """ returns Radar for successful parse, raises exception otherwise"""

    # if ... fields don't match
    #   raise ValueError("radar lines don't match")

    return Radar(
        # line1
        num=line1[:5].strip(),
        desig=line1[5:12].strip(),
        note1=line1[13].strip(),
        jdutc=date2JD(line1[15:32]),
        delay=float(line1[32:43] + '.' + line1[43:47]) if line1[42] != ' '
        else None,
        doppler=floatPx(line1[47:58] + '.' + line1[58:62]) if line1[57] != ' '
        else None,
        freq=float(line1[62:67] + '.' + line1[67:68] + line2[62:68]),
        txcod=line1[68:71],
        rxcod=line1[77:80],
        # line2
        ret=line2[32],
        udelay=float(line2[33:43] + '.' + line2[43:47]) if line2[42] != ' '
        else None,
        udoppler=float(line2[47:58] + '.' + line2[58:62]) if line2[57] != ' '
        else None
    )

def date2JD(dateObs):
    """
    Converts the date string to a JD floating point number, using the NOVAS routine
    The time scale of the returned value will be same as that of the input date.
    
    Parameters
    ----------
    RA      : string
    
    Returns
    -------
    degRA   : float
    
    Examples
    --------
    >>> ...
    
    """
    yr, mn, dy = parseDate(dateObs)
    if ' ' in dy:
        dy=dy.split(' ')[0]
    hr, dy = math.modf(float(dy))
    jd = novas.julian_date(int(yr), int(mn), int(dy), 24.*hr)
    return jd

def parseDate(dateObs):
    """
    Parses the date string from the 80-character record
    Called by "date2JD()"
    
    Parameters
    ----------
    dateObs      : string
    
    Returns
    -------
    yr, mn, dy   : strings
    
    Examples
    --------
    >>> ...
    
    """
    yr = dateObs[0:4]
    mn = dateObs[5:7]
    dy = dateObs[8:]
    return yr, mn, dy

def RA2hrRA(RA):
    hr = float(RA[0:2])
    if RA[5] == '.':
        return hr + float(RA[3:])/60. # decimal point in minutes
    mn = float(RA[3:5])
    try:
        sc = float(RA[6:])
    except ValueError:
        sc = 0
    hrRA = hr + 1./60. * (mn + 1./60. * sc)
    return hrRA

def Dec2degDec(Dec):
    """
    Convert the Dec string to float.
    
    Parameters
    ----------
    Dec      : string
    
    Returns
    -------
    degDec   : float
    units are in degrees
    
    Examples
    --------
    >>> ...
    
    """
    s = Dec[0]
    dg = float(Dec[1:3])
    if Dec[6] == '.':
        decDeg = dg + float(Dec[4:]) # decimal point in minutes
    else:
        mn = float(Dec[4:6])
        try:
            sc = float(Dec[7:])
        except ValueError:
            sc = 0
        degDec = dg + 1./60. * (mn + 1./60. * sc)
    if s == '-':
        degDec = -degDec
    return degDec


def RA2degRA(RA):
    """
    Convert the RA string to float.
    
    Parameters
    ----------
    RA      : string
    
    Returns
    -------
    degRA   : float
    units are in degrees
    
    Examples
    --------
    >>> ...
    
    """
    hr = float(RA[0:2])
    mn = float(RA[3:5])
    sc = float(RA[6:])
    degRA = 15.0*(hr + 1./60. * (mn + 1./60. * sc))
    return degRA


def floatMag(s):
    if s in ["     ", "  .  "]:
        return None
    return float(s)


def floatPx(s):
    f = float(s[1:])
    if s[0] == '-':
        return -f
    return f


def note2(line):
    return line[14]


def parse80(s):
    """generator, parses a text stream of 80 column observations.

    The parser is intended for relatively clean 80 column data.  It will
    quietly skip blank lines but otherwise does no reformatting before parsing.

    Parameters
    ----------
    s : io.TextIOBase
        text stream containing 80 column observations

    Yields
    ------
    object
        The generator yields an object for each observation or unparsable
        line.  Observations will be one of the named tuples defined above,
        Optical, SpaceBased, Radar, or Roving.  For unparsable lines the
        generator yields a tuple containing an Exception object and the
        unparseable line.
    """
    line1 = None
    for line in s:
        if line.isspace():
            continue
        try:
            if len(line) < 80:
                raise ValueError('< 80 characters')
            n2 = note2(line)
            if line1:
                l1 = line1
                line1 = None
                if ord(n2) != ord(note2(l1)) + 32:
                    raise ValueError('got {} want {}'.format(n2,
                        chr(ord(note2(l1)) + 32)))
                if n2 == 's':
                    yield parseSat(l1, line)
                elif n2 == 'r':
                    yield parseRad(l1, line)
                else:
                    yield parseRov(l1, line)
            elif n2 in 'SRV':
                line1 = line
            else:
                yield parseOpt(line)
        except Exception as e:
            yield e, line1, line


def run_file(txt_fn):
    """parse a file, write successful results on stdout, unsuccessful
    parse results to stderr
    """
    with open(txt_fn) as f:
        for o in parse80(f):
            if not isinstance(o[0], Exception):
                print(o)
            elif o[1] is None:
                print('obs80: {}: {}'.format(o[0], o[2].rstrip()),
                      file=sys.stderr)
            else:
                print("obs80: {} line1: {} line2: {}".format(*o),
                    file=sys.stderr)




# ---------------------------------------------------------------------------------------------------
# The subsequent Heliocentric & HCAnnotator were initially
# in a module named "obs80hc" by Sonia Keys.
# This was relocated here by Payne: Nov 2018
# There is an extremely heavy overlap between this code and
# Holman's "Observatory" class in "mpcutilities/mpc_library.
# This should be rationalized
# ---------------------------------------------------------------------------------------------------

# fields copied from Optical/SpaceBased + heliocentric position
# of observing site or spacecraft
Heliocentric = collections.namedtuple('Heliocentric',
                                      'num desig disc note1 note2 jdutc ra dec mag band cod hx hy yz')
"""named tuple with fields in common with the obs80 named tuple types Optical
    and SpaceBased, and annotated with heliocentric observing point coordinates.
    
    Parameters
    ----------
    num : string
        packed perm desig + comet/sat stuff
    desig : string
        packed provisional or temporary desig
    disc : string
        discovery asterisk
    note1 : string
        overloaded note/program code
    note2 : string
        overloaded record type/measurement technology
    jdutc : float
        time of obs parsed to JD
    ra : float
        RA parsed to hours
    dec : float
        dec parsed to degrees
    mag : float
        mag parsed to float
    band : string
        mag band
    cod : string
        3 character obs code
    hx : float
        heliocentric x in AU
    hy : float
        heliocentric y in AU
    hz : float
        heliocentric z in AU
"""




class HCAnnotator:
    """
        
    Annotates observations with heliocentric positions.
    
    [[ Payne: ]]
    
    [[ This was written by Sonia ]]
    
    [[ There is an extreme overlap with Holman's "Observatory" class in "mpcutilities/mpc_library]]
    
    [[ Needs to be rationalized ]]
    
    
    A JPL binary ephemeris file is required at the time the constructor is
    called.  There are apparently three ways to provide this file.
    
    1. Have the package novas_de405 installed and importable.
    2. Know the location of the binary ephemeris file and supply it as
    ephem_name to the constructor
    3. Know the location of the binary ephemeris file and set the environment
    variable EPHEMERIS_FILE.
    
    See readme, novas py readme and doc for novas eph_manager.ephem_open.
    
    Parameters
    ----------
    siteXYZ : dictionary
    dict as returned by module obscode, function siteXYZ.
    leapsec : object
    LeapSeconds object such as defined in module leapsec.
    ephem_name : string, optional
    File name of binary JPL ephemeris.
    """
    def __init__(self, siteXYZ, leapsec, ephem_name=None):
        eph_manager.ephem_open(ephem_name)
        self.siteXYZ = siteXYZ
        self.leapsec = leapsec
        self.cache = {}
    
    def annotate(self, parsed, other=None):
        """
        generator, annotates observations with heliocentric coordinates.
            
        Parameters
        ----------
        parsed : iterable
            iterable over objects including type Optical and SpaceBased,
            such as the generator obs80.read80
        other : list
            If a list is passed, objects other than Optical and SpaceBased
            are appended.  
            Pass None if the list is not of interest.
        
        Yields
        ------
        namedtuple
            Heliocentric named tuple for each Optical and SpaceBased
            observation, annotated with heliocentric positions.

        Example
        ------
        [[Payne:    I assume from Sonia's original ``if __name__ == '__main__':...'' ]]
        [[          that that the following usage is reasonable ... ]]
        >>> hc = HCAnnotator(
                            obscode.siteXYZ(obscode.read5c('obscode.dat')),
                            leapsec.LeapSeconds('leap-seconds.list'))
        >>> other = []
        >>> flat = hc.annotate(obs80.parse80(open(sys.argv[1])), other)
        >>> print(pandas.DataFrame(flat))
        >>> print(len(other), "others")
        >>> print("cache len:", len(hc.cache))

        """
        for obs in parsed:
            t = type(obs).__name__
            if t == 'Optical':
                hx, hy, hz = self._xyzOpt(obs)
                n2 = obs.note2
                cod = obs.cod
            elif t == 'SpaceBased':
                hx, hy, hz = self._xyzSat(obs)
                n2 = 'S'
                cod = obs.cod
            elif t == 'Roving':
                hx, hy, hz = self._xyzRov(obs)
                n2 = 'V'
                cod = '247'
            else:
                try:
                    other.append(obs)
                except AttributeError:
                    pass
                continue
            yield Heliocentric(obs.num, obs.desig, obs.disc,
                               obs.note1, n2, obs.jdutc, obs.ra, obs.dec,
                               obs.mag, obs.band, cod, hx, hy, hz)

    def _xyzOpt(self, obs):
        """
        calculates the heliocentric position of the observatory at UTC JD
        in equatorial cartesian coordinates.

        Note that for now it ignores polar motion, differences between TT and
        TDB, and differences between UTC and UT1.  The last is probably the
        biggest source of error.
        
        Payne: Compare / integrate with MPCU.Observatory.getObservatoryPosition()
        

        """
        key = (obs.cod, obs.jdutc)
        try:
            return self.cache[key]
        except KeyError:
            pass
        obsVec = self.siteXYZ[obs.cod]
        leaps = self.leapsec.getLeapSeconds(obs.jdutc)
        jd_tt = obs.jdutc + (32.184 + leaps) / (24.0 * 60 * 60)
        DUT1 = 0.0  # for now
        delta_t = 32.184 + leaps - DUT1
        jd_ut1 = obs.jdutc + DUT1
        xp, yp = 0.0, 0.0  # for now
        geocentric_vec = novas.ter2cel(jd_ut1, 0.0, delta_t, xp, yp, obsVec)
        geocentric_vec = PHYS.Rearth_AU * numpy.array(geocentric_vec)
        pos = MPCU.getEarthPosition(jd_tt)  # Might need to switch to tdb
        hc = (pos[0] + geocentric_vec[0],
              pos[1] + geocentric_vec[1],
              pos[2] + geocentric_vec[2])
        self.cache[key] = hc
        return hc

    def _xyzSat(self, obs):
        """ 
        calculates the heliocentric position of an observing spacecraft at
        UTC JD in equatorial cartesian coordinates.
        
        Payne: Compare / integrate with MPCU.Observatory.getSatellitePosition()
        """
        key = (obs.cod, obs.jdutc)
        try:
            return self.cache[key]
        except KeyError:
            pass
        leaps = self.leapsec.getLeapSeconds(obs.jdutc)
        jd_tt = obs.jdutc + (32.184 + leaps) / (24.0 * 60 * 60)
        pos = MPCU.getEarthPosition(jd_tt)  # Might need to switch to tdb
        hc = pos[0] + obs.x, pos[1] + obs.y, pos[2] + obs.z
        self.cache[key] = hc
        return hc

    def _xyzRov(self, obs):
        """ 
        calculates the heliocentric position of a roving observer at
        UTC JD in equatorial cartesian coordinates.
        
        Payne: This does *NOT* have an equivalent in MPCU.Observatory.
        
        """
        key = ('247', obs.jdutc)
        try:
            return self.cache[key]
        except KeyError:
            pass
        leaps = self.leapsec.getLeapSeconds(obs.jdutc)
        delta_t = 32.184 + leaps
        jd_tt = obs.jdutc + delta_t / (24.0 * 60 * 60)
        e = MPCU.getEarthPosition(jd_tt)  # Might need to switch to tdb
        observer = novas.make_observer_on_surface(obs.lat, obs.lon, obs.alt,
                  5, 990)
        # posvel returns GCRS
        o, _ = novas.geo_posvel(jd_tt, delta_t, observer)
        # add ICRS and GCRS?
        hc = (e[0] + o[0], e[1] + o[1], e[2] + o[2])
        self.cache[key] = hc
        return hc







if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("usage: obs80 <obs file>")
    # parse argument text file
    run_file(sys.argv[1])
