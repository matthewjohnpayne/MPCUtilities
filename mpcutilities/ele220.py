#!/usr/bin/env python3
"""

--------------------------------------------------------------

Nov 2018

- Payne

Originally written by Sonia Keys

Parsing and formatting for four orbit formats:

* the 220 column format for minor planets
* the 90 column format for minor planets
* the 220 column format for natural satellites of outer planets
* the 255 column format for comets

Parsing is provided with classes Ele220, Ele90, EleSat, EleComet.  

Formatting is provided class EleFields, but currently only implemented for the 220 and 90
column minor planet formats, not for the natural satellite or comet formats.

--------------------------------------------------------------

"""

import datetime
import math
import novas.compat as novas

__all__ = ['Ele220', 'Ele90', 'EleSat', 'EleComet']

k = 0.01720209895


_ps = {
    'd': 'DE200',
    'f': 'DE245',
    'h': 'DE403',
    'j': 'DE405'
}


class _Ele(str):
    """
        
    Base class for 220 and 90 character element lines.

    Parameters
    ----------
    str : string
        The 220 or character element line.  A trailing newline is okay.
        
    """

    # slice objects hold column bounds for the raw fields within the 220
    # character line.  Primarily for internal use, not as part of the API.
    _timePeri = slice(12, 24)
    _argPeri = slice(24, 34)
    _node = slice(34, 44)
    _inc = slice(44, 54)
    _periDist = slice(54, 64)
    _ecc = slice(64, 74)

    def timePeri(self):
        """
        Returns
        -------
        float
            Date of perihelion passage/JDT
        """
        return JDFromPacked(self[_Ele._timePeri])

    def argPeri(self):
        """
        Returns
        -------
        float
            Argument of perihelion (deg/J2000.0)
        """
        return float7(self[_Ele._argPeri])

    def node(self):
        """
        Returns
        -------
        float
            Longitude of ascending node (deg/J2000.0)
        """
        return float7(self[_Ele._node])

    def inc(self):
        """
        Returns
        -------
        float
            Inclination of orbit (deg/J2000.0)
        """
        return float7(self[_Ele._inc])

    def periDist(self):
        """
        Returns
        -------
        float
            Perihelion distance (AU)
        """
        return float9p(self[_Ele._periDist])

    def ecc(self):
        """
        Returns
        -------
        float
            Eccentricity
        """
        return float9(self[_Ele._ecc])


class Ele220(_Ele):
    """
        
    Holds a 220 character element line.

    Class Ele220 instantiates on a 220 character line, then methods interpret
    individual field values.  A design assumption is that many applications
    will not need all fields but only some subset.  Instantiation therefore
    does not do any parsing.  Individual methods interpret and return single
    field values.

    Thus depending on the application, Ele220 may not be best for orbit
    represention but only for parsing.  For example you might define an orbit
    class in some way suitable for computation, then instantiate an Ele220 on
    a line from an orbit file, then use the methods of the Ele220 object in
    constrution of your computation-suitable orbit object.  This way only the
    fields you need are parsed; they are parsed once and then the Ele220 object
    is no longer needed.

    Note the methods ma1, argPeri1, node1, inc1, ecc3, and a3 have little
    computational use.  They interpret fields appearing in the 220 column
    format as "readable" values with reduced precision.  Full precision
    values are either available as other methods or can be computed from
    other methods.

    Parameters
    ----------
    str : string
        220 character element line.  Extra length such as a trailing newline
        is okay.
    """
    _desig = slice(0, 7)
    _g = slice(7, 12)
    # 12:74 are standard elements
    _pertEpoch = slice(75, 80)
    _designation = slice(80, 90)
    _h = slice(91, 96)
    _maEpoch = slice(97, 102)
    _ma1 = slice(103, 108)
    _argPeri1 = slice(109, 114)
    _node1 = slice(115, 120)
    _inc1 = slice(121, 126)
    _ecc3 = slice(127, 132)
    _a3 = slice(133, 140)
    _numOpp = slice(141, 144)
    _arc = slice(144, 150)
    _numObs = slice(150, 156)
    _uNote = 157  # overloaded column
    _comp = slice(159, 168)
    _ref = slice(169, 178)
    _first = slice(179, 184)
    _last = slice(185, 190)
    _rms = slice(191, 196)
    _pertScheme = 197
    _pertCoarse = slice(199, 202)
    _pert = slice(203, 207)
    # there can be some capital letters after _pert.  no method.
    _numScore = slice(212, 214)
    _curOppScore = 214

    def desig(self):
        """
        Returns
        -------
        string
            a short designation, no longer than 7 characters, whitespace
            trimmed.  This may be a temporary or packed designation, either
            permanent or provisional.

        See also
        --------
        designation : Unpacked designation
        """
        return self[Ele220._desig].strip()

    def g(self):
        """
        Returns
        -------
        float
            Slope parameter, G
        """
        return float(self[Ele220._g])

    def pertEpoch(self):
        """
        Returns
        -------
        None or float
            float is epoch jd
        """
        pe = self[Ele220._pertEpoch]
        if pe.isspace():
            return None
        return JDFromPacked(pe)

    def designation(self):
        """
        Returns
        -------
        string
            a non-packed designation, up to 10 characters, whitespace
            trimmed.

        See also
        --------
        desig : Packed designation
        """
        return self[Ele220._designation].strip()

    def h(self):
        """
        Returns
        -------
        float
            Absolute magnitude, H
        """
        return float(self[Ele220._h])

    def maEpoch(self):
        """
        Returns
        -------
        float
            Epoch JD
        """
        return JDFromPacked(self[Ele220._maEpoch])

    def ma1(self):
        """
        Returns
        -------
        float
            Mean anomaly in degrees, precision limited to one decimal place.
        """
        return float(self[Ele220._ma1])

    def argPeri1(self):
        """
        Returns
        -------
        float
            Argument of perihelion in degrees, precision limited to one
            decimal place.
        """
        return float(self[Ele220._argPeri1])

    def node1(self):
        """
        Returns
        -------
        float
            Longitude of ascending node in degrees, precision limited to one
            decimal place.
        """
        return float(self[Ele220._node1])

    def inc1(self):
        """
        Returns
        -------
        float
            Inclination in degrees, precision limited to one decimal place.
        """
        return float(self[Ele220._inc1])

    def ecc3(self):
        """
        Returns
        -------
        float
            Eccentricity, precision limited to three decimal places.
        """
        return float(self[Ele220._ecc3])

    def a3(self):
        """
        Returns
        -------
        float
            Semi-major axis in AU, precision limited to three decimal places.
        """
        return float(self[Ele220._a3])

    def numOpp(self):
        """
        Returns
        -------
        int
            Number of oppositions included in solution
        """
        return int(self[Ele220._numOpp])

    def arc(self):
        """
        Returns
        -------
        int
            Arc length of observations included in solution (/days)
        """
        return int(self[Ele220._arc])

    def numObs(self):
        """
        Returns
        -------
        int
            Number of observations included in solution
        """
        return int(self[Ele220._numObs])

    def u(self):
        """
        Returns
        -------
        int
            U value (a number from 0 to 9) if present.
            10 if not present.  (10 is not a valid U value.)
        """
        try:
            return int(self[Ele220._uNote])
        except ValueError:
            return 10

    def uNote(self):
        """
        Returns
        -------
        string
            single character note from the U field, if one exists, '' otherwise
        """
        n = self[Ele220._uNote]
        if n == ' ' or n.isdigit():
            return ''
        return n

    def comp(self):
        """
        Returns
        -------
        string
            Computer name, up to 9 characters, whitespace trimmed
        """
        return self[Ele220._comp].strip()

    def ref(self):
        """
        Returns
        -------
        string
            Reference, up to 9 characters, whitespace trimmed
        """
        return self[Ele220._ref].strip()

    def first(self):
        """
        Returns
        -------
        float
            Date of first included observation as JD
        """
        return JDFromPacked(self[Ele220._first])

    def last(self):
        """
        Returns
        -------
        float
            Date of last included observation as JD
        """
        return JDFromPacked(self[Ele220._last])

    def rms(self):
        """
        Returns
        -------
        float
            r.m.s. fit of included obserations (/")
        """
        return float(self[Ele220._rms])

    def pertScheme(self):
        """
        Returns
        -------
        string
            a string such as 'DE200', 'DE403', 'DE405', or ''
        """
        if self[Ele220._pertEpoch.start] == ' ':
            return ''
        return _ps.get(self[Ele220._pertScheme], '')

    def pertCoarse(self):
        """
        Returns
        -------
        string
            a string such as 'M-v' or ''
        """
        if self[Ele220._pertEpoch.start] == ' ':
            return ''
        return self[Ele220._pertCoarse].strip()

    def perturbers(self):
        """
        Returns
        -------
        string
            4 character perturber code if present, '' otherwise
        """
        if self[Ele220._pertEpoch.start] == ' ':
            return ''
        return self[Ele220._pert]

    def numScore(self):
        """
        Returns
        -------
        int
            a number from 0 to 219 if present.
            -1 if not present.
        """
        s = self[Ele220._numScore]
        if s.isspace():
            return -1
        n = int(s[1])
        return n if s[0] == ' ' else b62(s[0]) * 10 + n

    def curOppScore(self):
        """
        Returns
        -------
        int
            a number from 0 to 9 if present.
            -1 if not present.
        """
        s = self[Ele220._curOppScore]
        return -1 if s == ' ' else int(s)


class Ele90(_Ele):
    """holds a 90 character element line.

    Parameters
    ----------
    str : string
        The 90 character element line.  A trailing newline is okay.

    See also
    --------
    Ele220
    """

    _desig = slice(0, 7)
    _h = slice(7, 12)
    # 12:74 are standard elements
    _first = slice(80, 85)
    _last = slice(85, 90)

    def desig(self):
        """
        Returns
        -------
        string
            a short designation, no longer than 7 characters, whitespace
            trimmed.  This may be a temporary or packed designation, either
            permanent or provisional.

        See also
        --------
        designation : Unpacked designation
        """
        return self[Ele90._desig].strip()

    def h(self):
        """
        Returns
        -------
        float
            Absolute magnitude, H, if field is filled
        None
            if field is blank
        """
        try:
            return float(self[Ele90._h])
        except ValueError:
            if not self[Ele90._h].isspace():
                raise

    def first(self):
        """
        Returns
        -------
        float
            Date of first included observation as JD
        """
        return JDFromPacked(self[Ele90._first])

    def last(self):
        """
        Returns
        -------
        float
            Date of last included observation as JD
        """
        return JDFromPacked(self[Ele90._last])


class EleSat(_Ele):
    """
        
    Holds a 220 character element line for a natural satellite

    Class EleSat instantiates on a 220 character line, then methods interpret
    individual field values.  A design assumption is that many applications
    will not need all fields but only some subset.  Instantiation therefore
    does not do any parsing.  Individual methods interpret and return single
    field values.

    Parameters
    ----------
    str : string
        220 character element line.  Extra length such as a trailing newline
        is okay.
    """
    _permDesig = slice(0, 4)
    # [4] is always 'S'
    _tempDesig = slice(5, 12)
    # 12:74 are standard elements
    _pertEpoch = slice(75, 80)
    _designation = slice(80, 91)
    _h = slice(91, 96)
    _maEpoch = slice(97, 102)
    _ma1 = slice(103, 108)
    _argPeri1 = slice(109, 114)
    _node1 = slice(115, 120)
    _inc1 = slice(121, 126)
    _ecc3 = slice(127, 132)
    _a3 = slice(133, 140)
    _arc = slice(144, 150)
    _numObs = slice(150, 156)
    _comp = slice(159, 168)
    _ref = slice(169, 178)
    _first = slice(179, 184)
    _last = slice(185, 190)
    _rms = slice(191, 196)
    _pertScheme = 197
    _pertCoarse = slice(199, 202)
    _pert = slice(203, 207)
    _numScore = slice(212, 214)
    _curOppScore = 214

    def permDesig(self):
        """
        Returns
        -------
        string
            permanent designation, 4 characters if present, '' if not present.
        """
        return self[EleSat._permDesig].strip()

    def tempDesig(self):
        """
        Returns
        -------
        string
            temporary designation, 7 characters if present, '' if not present.
        """
        return self[EleSat._tempDesig].strip()

    def pertEpoch(self):
        """
        Returns
        -------
        float
            epoch jd
        """
        return JDFromPacked(self[EleSat._pertEpoch])

    def designation(self):
        """
        Returns
        -------
        string
            unpacked designation.
        """
        return self[EleSat._designation].strip()

    def h(self):
        """
        Returns
        -------
        float
            Absolute magnitude, H
        """
        return float(self[EleSat._h])

    def maEpoch(self):
        """
        Returns
        -------
        float
            Epoch JD
        """
        return JDFromPacked(self[EleSat._maEpoch])

    def ma1(self):
        """
        Returns
        -------
        float
            Mean anomaly in degrees, precision limited to one decimal place.
        """
        return float(self[EleSat._ma1])

    def argPeri1(self):
        """
        Returns
        -------
        float
            Argument of perihelion in degrees, precision limited to one
            decimal place.
        """
        return float(self[EleSat._argPeri1])

    def node1(self):
        """
        Returns
        -------
        float
            Longitude of ascending node in degrees, precision limited to one
            decimal place.
        """
        return float(self[EleSat._node1])

    def inc1(self):
        """
        Returns
        -------
        float
            Inclination in degrees, precision limited to one decimal place.
        """
        return float(self[EleSat._inc1])

    def ecc3(self):
        """
        Returns
        -------
        float
            Eccentricity, precision limited to three decimal places.
        """
        return float(self[EleSat._ecc3])

    def a3(self):
        """
        Returns
        -------
        float
            Semi-major axis in AU, precision limited to three decimal places.
        """
        return float(self[EleSat._a3])

    def arc(self):
        """
        Returns
        -------
        int
            Arc length of observations included in solution
        """
        return int(self[EleSat._arc])

    def numObs(self):
        """
        Returns
        -------
        int
            Number of observations included in solution
        """
        return int(self[EleSat._numObs])

    def comp(self):
        """
        Returns
        -------
        string
            Computer name, up to 9 characters, whitespace trimmed
        """
        return self[EleSat._comp].strip()

    def ref(self):
        """
        Returns
        -------
        string
            Reference, up to 9 characters, whitespace trimmed
        """
        return self[EleSat._ref].strip()

    def first(self):
        """
        Returns
        -------
        float
            Date of first included observation as JD
        """
        return JDFromPacked(self[EleSat._first])

    def last(self):
        """
        Returns
        -------
        float
            Date of last included observation as JD
        """
        return JDFromPacked(self[EleSat._last])

    def rms(self):
        """
        Returns
        -------
        float
            r.m.s. fit of included obserations
        """
        return float(self[EleSat._rms])

    def pertScheme(self):
        """
        Returns
        -------
        string
            a string such as 'DE200', 'DE403', 'DE405', or ''
        """
        return _ps.get(self[EleSat._pertScheme], '')

    def pertCoarse(self):
        """
        Returns
        -------
        string
            a string such as 'M-v' or ''
        """
        return self[EleSat._pertCoarse].strip()

    def perturbers(self):
        """
        Returns
        -------
        string
            4 character perturber code if present, '' otherwise
        """
        return self[EleSat._pert]

    def numScore(self):
        """
        Returns
        -------
        int
            a number from 0 to 219 if present.
            -1 if not present.
        """
        s = self[EleSat._numScore]
        if s.isspace():
            return -1
        n = int(s[1])
        return n if s[0] == ' ' else b62(s[0]) * 10 + n

    def curOppScore(self):
        """
        Returns
        -------
        int
            a number from 0 to 9 if present.
            -1 if not present.
        """
        s = self[EleSat._curOppScore]
        return -1 if s == ' ' else int(s)


class EleComet(_Ele):
    """holds a 256 character comet element line.

    Parameters
    ----------
    str : string
        The 256 character element line.  A trailing newline is okay.
    """

    _permDesig = slice(0, 5)   # includes orbit type at [4]
    _tempDesig = slice(4, 12)  # includes orbit type at [4]
    _tempDes1 = 5  # first position after the orbit type
    _frag = 11
    _frag2 = slice(10, 12)
    # 12:74 are standard elements
    # [74] is usually blank but sometimes has an 'a', or an 'S'.
    # there is no method for this.
    _pertEpoch = slice(75, 80)
    _name = slice(80, 109)
    _comp = slice(109, 120)  # there can be a leading @
    _pubYear = slice(121, 125)
    _pub = slice(125, 129)
    # unknown fields here.  looks like float, int, float, int, and a
    # date formatted as a string.
    _timeScale = 163
    _pubVol = slice(165, 174)
    _numObs = slice(176, 181)  # first char can be A or N, last can be +
    _nonGravStar = 181
    _first = slice(183, 188)
    _last = slice(189, 194)
    _pertScheme = 195
    # there can be a capital letter in 196.  no method.
    _pertCoarse = slice(197, 200)
    _pert = slice(201, 205)
    # [206] typically has a digit in it.
    # 208, 209, and 210 are note fields,  no methods for these, except
    # that a ! in 210 obviously means that fields in [183:210] cannot be
    # interpreted by _first, _last, _pertScheme, _pertCoarse, _pert.
    # if 210 contains a !, [183:210] has unknown meaning.
    _notes = slice(208, 210)
    # a '!' sometimes appears in [219]
    # [222:251] hold non-grav parameters, at least if [_nonGravStar] is '*'.
    # Interpretation of these columns in non-* cases is unknown.
    _ng_style = 221
    _ng_a1 = slice(222, 232)
    _ng_a2 = slice(232, 243)
    _ng_a3 = slice(243, 251)
    _rms = slice(251, 255)

    def permDesig(self):
        """
        Returns
        -------
        string
            permanent designation, 5 characters if present, '' if not.
        """
        return self[EleComet._permDesig] if self[0] > ' ' else ''

    def tempDesig(self):
        """
        Returns
        -------
        string
            temporary designation, including leading orbit type character,
            8 characters altogether if present, '' if not present.
        """
        if self[EleComet._tempDes1] > ' ':
            return self[EleComet._tempDesig]
        return ''

    def frag(self):
        """
        Returns
        -------
        string
            1 or 2 character fragment identifer if present.
            '' if not present.
        """
        # _frag is the last position of _tempDesig. if there's no lowercase
        # letter there, there's no fragment identifier.
        if self[EleComet._frag] < 'a':
            return ''
        # non-space after the orbit type character at the start of tempDesig
        # means the tempDesig is present and the fragment can only be the
        # single character.
        if self[EleComet._tempDesig.start + 1] > ' ':
            return self[EleComet._frag]
        # else, with no temp des, the field is free to hold two character
        # fragment designations.
        if self[EleComet._frag2.start] == ' ':
            return self[EleComet._frag]
        return self[EleComet._frag2]

    def pertEpoch(self):
        """
        Returns
        -------
        None or float
            float is epoch jd
        """
        pe = self[EleComet._pertEpoch]
        if pe.isspace():
            return None
        return JDFromPacked(pe)

    def name(self):
        """
        Returns
        -------
        string
            Comet name
        """
        return self[EleComet._name].strip()

    def comp(self):
        """
        Returns
        -------
        string
            Orbit computer name.  The first character can be '@', which
            means it can be expanded with an external list,
            comet_computers.txt.
        """
        return self[EleComet._comp].strip()

    def pubYear(self):
        """
        Returns
        -------
        string
        """
        return self[EleComet._pubYear].strip()

    def pub(self):
        """
        Returns
        -------
        string
        """
        return self[EleComet._pub].strip()

    def timeScale(self):
        """
        Returns
        -------
        string
            Either "TT/TDT" or "UTC".
        """
        if self[EleComet._timePeri.start] < 'J':
            return "TT/TDT" if self[EleComet._timeScale] == 'T' else "UTC"
        return "UTC" if self[EleComet._timeScale] == 'U' else "TT/TDT"

    def pubVol(self):
        """
        Returns
        -------
        string
        """
        return self[EleComet._pubVol].strip()

    def numObs(self):
        """
        Returns
        -------
        int
            Number of observations used in orbit solution.
        """
        n = self[EleComet._numObs]
        if n[0] in "AN":
            return int(n[1:])
        if n[4] == '+':
            return int(n[:4])
        return int(n)

    def numObsAN(self):
        """
        Returns
        -------
        string
            'A' or 'N' where one of these characters appears before numObs.
        """
        n0 = self[EleComet._numObs.start]
        return n0 if n0 in "AN" else ""

    def numObsPlus(self):
        """
        Returns
        -------
        bool
            True if a '+' appears following numObs.
        """
        return self[EleComet._numObs.stop - 1] == '+'

    def nonGrav(self):
        """
        Returns
        -------
        None or tuple
            None if no non-grav parameters present, otherwise:
            3-tuple of (style, a1, a2), or
            4-tuple of (style, a1, a2, a3)
            where style is a single character and a1-a3 are floats
        """
        if self[EleComet._nonGravStar] != '*':
            return None
        st = self[EleComet._ng_style]
        a1 = float(self[EleComet._ng_a1])
        a2 = float(self[EleComet._ng_a2])
        s3 = self[EleComet._ng_a3].strip()
        if s3 == '':
            return st, a1, a2
        return st, a1, a2, float(s3)

    def first(self):
        """
        Returns
        -------
        float or None
            Date of first included observation as JD if valid, None otherwise
        """
        if self[210] == '!':
            return None
        return JDFromPacked(self[EleComet._first].strip())

    def last(self):
        """
        Returns
        -------
        float or None
            Date of last included observation as JD if valid, None otherwise
        """
        if self[210] == '!':
            return None
        return JDFromPacked(self[EleComet._last].strip())

    def pertScheme(self):
        """
        Returns
        -------
        string
            a string such as 'DE200', 'DE403', 'DE405', or ''
        """
        if self[EleComet._pertEpoch.start] == ' ' or self[210] == '!':
            return ''
        return _ps.get(self[EleComet._pertScheme], '')

    def pertCoarse(self):
        """
        Returns
        -------
        string
            a string such as 'M-v' or ''
        """
        if self[EleComet._pertEpoch.start] == ' ' or self[210] == '!':
            return ''
        return self[EleComet._pertCoarse].strip()

    def perturbers(self):
        """
        Returns
        -------
        string
            4 character perturber code if present, '' otherwise
        """
        if self[EleComet._pertEpoch.start] == ' ' or self[210] == '!':
            return ''
        return self[EleComet._pert]

    def notes(self):
        """
        Returns
        -------
        string
            up to two note codes.
        """
        return self[EleComet._notes].strip()

    def rms(self):
        """
        Returns
        -------
        float
            r.m.s. fit of included obserations if rms is present.
            -1 if rms is not present.
        """
        s = self[EleComet._rms]
        return -1 if s.isspace() else float(s)


_i62 = [None] * (ord('z') + 1)
for i in range(10):
    _i62[ord('0') + i] = i
for i in range(26):
    _i62[ord('A') + i] = i + 10
for i in range(26):
    _i62[ord('a') + i] = i + 10 + 26


def b62(c):
    try:
        n = _i62[ord(c)]
        if n is not None:
            return n
    except IndexError:
        pass
    raise ValueError(c)


def unpackDate(s):
    """
    Returns
    -------
    int
        year
    int
        month
    int
        day
    float
        fraction of day
    """
    if len(s) < 3:  # allow as little as a packed year
        raise ValueError(s)

    # century
    c = ord(s[0])
    if c >= ord('A'):
        c = 10 + c - ord('A')
    else:
        c -= ord('0')

    if c > 25:  # allow a few centuries into the future,
        raise ValueError(s)  # call anything beyond that an error

    try:
        y = int(s[1:3])
    except ValueError:
        raise ValueError(s)

    if c >= 0:
        y = c * 100 + y
    else:
        y = (c + 1) * 100 - y

    m = 1
    if len(s) > 3:
        try:
            m = b62(s[3])
        except ValueError:
            raise ValueError(s)
        if m == 0 or m > 12:
            raise ValueError(s)

    d = 1
    if len(s) > 4:
        try:
            d = b62(s[4])
        except ValueError:
            raise ValueError(s)
        if d == 0 or d > 31:
            raise ValueError(s)

    frac = 0.
    if len(s) > 5:
        try:
            frac = float('.' + s[5:])
        except ValueError:
            raise ValueError(s)

    return y, m, d, frac


def JDFromPacked(s):
    y, m, d, frac = unpackDate(s)
    return to_julian_date(y, m, d + frac)


def to_julian_date(year, month, day):
    frac, d = math.modf(day)
    j = novas.julian_date(year, month, int(d), frac * 24)
    if j >= 2299159.5:  # start of Gregorian calendar, October 15, 1582
        return j
    # else compute correction for Julian calendar
    y = year - (14 - month) // 12
    return j + y // 100 - y // 400 - 2


def packDate(j, w):
    y, m, d = from_julian_date(j)
    prec = w - 5
    # split fraction from day
    frac, d = math.modf(d)
    d = int(d)
    # round fraction to target precision
    p = 10. ** prec
    frac = math.floor(p * frac + .5)
    if frac == p:  # round up to next day
        d = datetime.date(y, m, d) + datetime.timedelta(days=1)
        y = d.year
        m = d.month
        d = d.day
        frac = 0
    # after rounding check year can be packed
    if y < 1800 or y > 2099:
        raise ValueError('cannot pack date {}'.format(j))
    if prec == 0:
        frac = ''
    else:
        frac = '%0*d' % (prec, frac)
    if d < 10:
        d = str(d)
    else:
        d = chr(int(d) - 10 + ord('A'))
    return '{}{:02d}{:X}{}{}'.format(chr(y // 100 - 18 + ord('I')), y % 100,
        m, d, frac)


def from_julian_date(j):
    if j < 2299159.5:  # start of Gregorian calendar, October 15, 1582
        # this could be implemented if needed.  not sure it is needed.
        # one possible implementation would be to compute a correction as
        # was done for to_julian_date.
        raise ValueError('JD before start of Gregorian calendar unsupported')
    y, m, d, hr = novas.cal_date(j)
    return y, m, d + hr / 24


def float7(s):
    return int(s) * 1e-7


def float9(s):
    return int(s) * 1e-9


def float9p(s):
    if '.' in s:
        return float(s)
    return int(s) * 1e-9


class EleFields:
    """
    
    Abstract base class providing formatting methods.

    Formatting methods fmt220 and fmt90 format EleFields attributes but
    the attributes are not part of the EleFields class definition.  
    It may be most convenient to use EleFields as an abstract base class 
    for a class that defines the needed attributes.  
    Alternatively the constructor can define attributes.

    Parameters
    ----------
    kw : keywords
        keyword parameters will be set as attributes.
    """
    def __init__(self, **kw):
        self.__dict__.update(kw)

    def fmt220(self):
        """
        formats attributes into a string in the 220 character format.

        Recognized attributes:

        desig, g, timePeri, argPeri, node, inc, periDist, ecc, designation,
        h, maEpoch, ma, a, numOpp, arc, numObs, comp, ref, first, last, rms.

        Fields of the 220 character format without defined attributes are
        left blank.  Note there are not generally separate high and low
        precision attributes.  Mean anomaly ma and semi-major axis a, only
        being formatted in low precision, may be derived from high precision
        attributes such as periDist and ecc.  The reverse computations are
        not automatically performed however.  No assumption is made on the
        precision of attributes ma and a and so no high precision fields will
        be derived from them.

        Returns
        -------
        string
            formatted 220 character string

        Examples
        --------
        >>> # use EleFields as an abstract base class
        >>> class aei(EleFields):
        ...   def __init__(self, a, e, i):
        ...     self.a = a
        ...     self.ecc = e
        ...     self.inc = i
        ...
        >>> x = aei(3.1, .1, 20)
        >>> line = x.fmt220()
        >>> line[Ele220._inc]
        '0200000000'
        >>> # periDist not automatically computed from a and ecc.
        >>> line[Ele220._periDist]
        '          '
        >>> line[Ele220._ecc]
        '0100000000'
        >>> line
        '                                            0200000000          0100000000                                                20.0 0.100   3.100                                                                                '

        >>> # construct with kwargs
        >>> line = EleFields(periDist=2.79, ecc=.1, inc=20).fmt220()
        >>> line[Ele220._periDist]
        '2790000000'
        >>> line[Ele220._a3]  # a3 computed from periDist and ecc
        '  3.100'
        >>> line
        '                                            020000000027900000000100000000                                                20.0 0.100   3.100                                                                                '
        """
        return ''.join([
            self._desig(),
            self._g(),
            self._timePeri(),
            self._argPeri(),
            self._node(),
            self._inc(),
            self._periDist(),
            self._ecc(),
            '      ',  # undocumented field, appears to be a packed date.
            self._designation(),
            ' ',
            self._h(),
            '  ',
            self._maEpoch(),
            ' ',
            self._ma1(),
            ' ',
            self._argPeri1(),
            ' ',
            self._node1(),
            ' ',
            self._inc1(),
            ' ',
            self._ecc3(),
            ' ',
            self._a3(),
            ' ',
            self._numOpp(),
            ' ',
            self._arc(),
            ' ',
            self._numObs(),
            '   ',  # undocumented field
            self._comp(),
            ' ',
            self._ref(),
            ' ',
            self._first(),
            ' ',
            self._last(),
            ' ',
            self._rms(),
            '                        '])  # multiple undocumented fields

    def _desig(self):
        # all field formatters roughly follow this pattern:
        try:
            v = self.desig  # first just attempt to access the field,
        except AttributeError:
            return '       '  # and return blank if it's not there.
        if v is None:
            return '       '  # also return blank if it's None.
        if len(v) <= 7:  # perform any range check
            return v.ljust(7)  # return the formatted field if it's okay
        raise ValueError('desig > 7 characters')  # raise ValueError otherwise

    def _g(self):
        try:
            v = self.g
        except AttributeError:
            return '     '
        if v is None:
            return '     '
        # 0 to .4 is the typical range for g.  the limits here hopefully
        # are enough to allow some unusual yet intentional values.
        if v > -.5 and v < 1:
            return '{:5.2f}'.format(v)
        raise ValueError('invalid g')

    def _timePeri(self):
        try:
            v = self.timePeri
        except AttributeError:
            return '            '
        if v is None:
            return '            '
        try:
            return packDate(v, 12)
        except ValueError:
            raise ValueError('invalid timePeri')

    def _argPeri(self):
        try:
            v = self.argPeri
        except AttributeError:
            return '          '
        if v is None:
            return '          '
        s = '{:011.7f}'.format(v % 360.)
        if s == '360.0000000':
            return '0000000000'
        return s[:3] + s[4:]

    def _node(self):
        try:
            v = self.node
        except AttributeError:
            return '          '
        if v is None:
            return '          '
        s = '{:011.7f}'.format(v % 360.)
        if s == '360.0000000':
            return '0000000000'
        return s[:3] + s[4:]

    def _inc(self):
        try:
            v = self.inc
        except AttributeError:
            return '          '
        if v is None:
            return '          '
        if v >= 0 and v <= 180:
            s = '{:011.7f}'.format(v)
            return s[:3] + s[4:]
        raise ValueError('invalid inc')

    def _periDist(self):
        try:
            v = self.periDist
        except AttributeError:
            return '          '
        if v is None:
            return '          '
        if v > 0:
            s = '{:011.9f}'.format(v)
            if len(s) == 11:
                return s[:1] + s[2:]
            return '%.*f' % (19 - len(s), v)
        raise ValueError('invalid periDist')

    def _ecc(self):
        try:
            v = self.ecc
        except AttributeError:
            return '          '
        if v is None:
            return '          '
        if v >= 0 and v < 1.1:
            s = '{:011.9f}'.format(v)
            return s[:1] + s[2:]
        raise ValueError('invalid ecc')

    def _designation(self):
        try:
            v = self.designation
        except AttributeError:
            return '          '
        if v is None:
            return '          '
        if len(v) <= 10:
            return v.ljust(10)
        raise ValueError('designation > 10 characters')

    def _h(self):
        try:
            v = self.h
        except AttributeError:
            return '    '
        if v is None:
            return '    '
        if v > -1 and v < 60:
            return '{:4.1f}'.format(v)
        raise ValueError('invalid h')

    def _maEpoch(self):
        try:
            v = self.maEpoch
        except AttributeError:
            return '     '
        if v is None:
            return '     '
        try:
            return packDate(v, 5)
        except ValueError:
            raise ValueError('invalid maEpoch')

    def _ma1(self):
        try:
            v = self.ma
        except AttributeError:
            try:
                a = self.periDist / (1 - self.ecc)
                n = k / math.sqrt(a * a * a) * 180 / math.pi
                v = n * (self.maEpoch - self.timePeri)
            except AttributeError:
                return '     '
        if v is None:
            return '     '
        s = '{:5.1f}'.format(v % 360)
        if s == '360.0':
            return '  0.0'
        return s

    def _argPeri1(self):
        try:
            v = self.argPeri
        except AttributeError:
            return '     '
        if v is None:
            return '     '
        s = '{:5.1f}'.format(v % 360.)
        if s == '360.0':
            return '  0.0'
        return s

    def _node1(self):
        try:
            v = self.node
        except AttributeError:
            return '     '
        if v is None:
            return '     '
        s = '{:5.1f}'.format(v % 360.)
        if s == '360.0':
            return '  0.0'
        return s

    def _inc1(self):
        try:
            v = self.inc
        except AttributeError:
            return '     '
        if v is None:
            return '     '
        if v >= 0 and v <= 180:
            return '{:5.1f}'.format(v)
        raise ValueError('invalid inc')

    def _ecc3(self):
        try:
            v = self.ecc
        except AttributeError:
            return '     '
        if v is None:
            return '     '
        if v >= 0 and v < 1.1:
            return '{:5.3f}'.format(v)
        raise ValueError('invalid ecc')

    def _a3(self):
        try:
            v = self.a
        except AttributeError:
            try:
                v = self.periDist / (1 - self.ecc)
            except AttributeError:
                return '       '
        if v is None:
            return '       '
        if v > 0 and v <= 999.999:
            return '{:7.3f}'.format(v)
        # if a doesn't fit in the format, just leave it blank.
        return '       '

    def _numOpp(self):
        try:
            v = self.numOpp
        except AttributeError:
            return '   '
        if v is None:
            return '   '
        # numOpp must be strictly > 0.  300 is still a quite generous upper
        # limit while disallowing any use of 999 to mean anything special.
        if v > 0 and v < 300:
            return '{:3d}'.format(v)
        raise ValueError('invalid numOpp')

    def _arc(self):
        try:
            v = self.arc
        except AttributeError:
            return '     '
        if v is None:
            return '     '
        # arc can be 0 in the case of rounding an arc < 12 hours.
        # 30000 a generous upper limit disallowing 99999.
        if v >= 0 and v <= 30000:
            return '{:5d}'.format(v)
        raise ValueError('invalid arc')

    def _numObs(self):
        try:
            v = self.numObs
        except AttributeError:
            return '     '
        if v is None:
            return '     '
        # numObs must be strictly > 0.
        # 30000 a generous upper limit disallowing 99999.
        if v > 0 and v < 30000:
            return '{:5d}'.format(v)
        raise ValueError('invalid numObs')

    def _comp(self):
        try:
            v = self.comp
        except AttributeError:
            return '         '
        if v is None:
            return '         '
        if len(v) <= 9:
            return v.ljust(9)
        raise ValueError('comp > 9 characters')

    def _ref(self):
        try:
            v = self.ref
        except AttributeError:
            return '         '
        if v is None:
            return '         '
        if len(v) <= 9:
            return v.ljust(9)
        raise ValueError('ref > 9 characters')

    def _first(self):
        try:
            v = self.first
        except AttributeError:
            return '     '
        if v is None:
            return '     '
        try:
            return packDate(v, 5)
        except ValueError:
            raise ValueError('invalid first')

    def _last(self):
        try:
            v = self.last
        except AttributeError:
            return '     '
        if v is None:
            return '     '
        try:
            return packDate(v, 5)
        except ValueError:
            raise ValueError('invalid last')

    def _rms(self):
        try:
            v = self.rms
        except AttributeError:
            return '     '
        if v is None:
            return '     '
        if v >= 0 and v < 60:
            return '{:5.2f}'.format(v)
        raise ValueError('invalid rms')

    def fmt90(self):
        """
        formats attributes into a string in the 90 character format.

        Recognized attributes:

        desig, h, timePeri, argPeri, node, inc, periDist, ecc, first, last.

        Returns
        -------
        string
            formatted 90 character string
        """
        return ''.join([
            self._desig(),
            self._h(),
            ' ',
            self._timePeri(),
            self._argPeri(),
            self._node(),
            self._inc(),
            self._periDist(),
            self._ecc(),
            '      ',  # unknown purpose
            self._first(),
            self._last()])
