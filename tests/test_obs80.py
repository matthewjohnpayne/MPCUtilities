# mpcutilities/tests/test_obs80.py

"""

--------------------------------------------------------------

Oct 2018

Payne copied in tests primarily written by Keys

Test the obs80 parse functions


--------------------------------------------------------------

"""
# Import third-party packages
# --------------------------------------------------------------
import numpy
import astropy.coordinates

# Importing of local modules/packages required for this test
# --------------------------------------------------------------
#import mpcutilities.phys_const  as PHYS
#import mpcutilities.classes     as Classes


# Import the specific package/module/function we are testing
# --------------------------------------------------------------
import mpcutilities.obs80 as obs80


# Basic test(s)
# ----------------------------------------------------------

def test_o1():
    o = obs80.parseOpt('00433         A1893 10 29.4132  06 08 59.32 +53 39 04.2                 HA053802')  # noqa E501
    assert o.num == '00433'
    assert o.desig == ''
    assert o.disc == ''
    assert o.note1 == ''
    assert o.note2 == 'A'
    t = astropy.time.Time(o.jdutc, format='jd').datetime
    assert t.year == 1893
    assert t.month == 10
    assert t.day == 29
    assert numpy.isclose(.4132, (t.hour + (t.minute + (t.second +
        t.microsecond / 1e6) / 60) / 60) / 24)
    h, m, s = astropy.coordinates.Angle(o.ra, unit='h').hms
    assert h == 6
    assert m == 8
    assert numpy.isclose(s, 59.32)
    sign, d, m, s = astropy.coordinates.Angle(o.dec, unit='deg').signed_dms
    assert sign == 1
    assert d == 53
    assert m == 39
    assert numpy.isclose(s, 4.2)
    assert o.mag is None
    assert o.band == ''
    assert o.cod == "802"


def test_o2():
    """negative dec"""
    o = obs80.parseOpt('00433         A1894 02 17.1174  07 36 21.35 -00 38 13.2                 HA053802')  # noqa E501
    sign, d, m, s = astropy.coordinates.Angle(o.dec, unit='deg').signed_dms
    assert sign == -1
    assert d == 0
    assert m == 38
    assert numpy.isclose(s, 13.2)


def test_o3():
    """decimal minutes in RA, missing seconds in dec, mag, no band"""
    o = obs80.parseOpt('00433         A1898 08 13.89338 21 33.8     -05 58               10.0   AN147020')  # noqa E501
    h, m, s = astropy.coordinates.Angle(o.ra, unit='h').hms
    assert h == 21
    assert numpy.isclose(m + s / 60., 33.8)
    sign, d, m, s = astropy.coordinates.Angle(o.dec, unit='deg').signed_dms
    assert sign == -1
    assert d == 5
    assert m == 58
    assert numpy.isclose(s, 0)
    assert o.mag == 10
    assert o.band == ''


def test_o4():
    """one decimal place in sec of ra, no decimal places in sec dec"""
    o = obs80.parseOpt('00433         A1898 08 14.96732 21 31 50.8  -05 57 43            11.0   AN147537')  # noqa E501
    h, m, s = astropy.coordinates.Angle(o.ra, unit='h').hms
    assert h == 21
    assert m == 31
    assert numpy.isclose(s, 50.8)
    sign, d, m, s = astropy.coordinates.Angle(o.dec, unit='deg').signed_dms
    assert sign == -1
    assert d == 5
    assert m == 57
    assert numpy.isclose(s, 43)


def test_o5():
    """empty note2"""
    o = obs80.parseOpt('00433          1898 09 20.88897 20 42 53.05 -05 59 42.4                m24429000')  # noqa E501
    assert o.note2 == ''


def test_o6():
    """both perm and prov, disc *, very coarse time, ra, dec, one dp in mag"""
    o = obs80.parseOpt('00433J56P00C* A1956 08 07.14    19 43.6     -20 47               14.4   M1581760')  # noqa E501
    assert o.num == '00433'
    assert o.desig == 'J56P00C'
    assert o.disc == '*'
    t = astropy.time.Time(o.jdutc, format='jd').datetime
    assert t.year == 1956
    assert t.month == 8
    assert t.day == 7
    assert numpy.isclose(.14, (t.hour + (t.minute + (t.second +
        t.microsecond / 1e6) / 60) / 60) / 24)
    h, m, s = astropy.coordinates.Angle(o.ra, unit='h').hms
    assert h == 19
    assert numpy.isclose(m + s / 60., 43.6)
    sign, d, m, s = astropy.coordinates.Angle(o.dec, unit='deg').signed_dms
    assert sign == -1
    assert d == 20
    assert m == 47
    assert numpy.isclose(s, 0)
    assert o.mag == 14.4


def test_o7():
    """note2"""
    o = obs80.parseOpt('00433         E1975 01 24.01889007 44 23.906+24 23 59.87                v4566244')  # noqa E501
    assert o.note2 == 'E'


def test_o8():
    """note1"""
    o = obs80.parseOpt('00433        4A1979 08 27.24896 21 15 51.51 -05 15 55.5                 M5026688')  # noqa E501
    assert o.note1 == '4'


def test_o9():
    """note2, extra precision, band"""
    o = obs80.parseOpt('00433         T1988 10 05.04698200 52 27.557+38 16 05.91         10.74V CMC05950')  # noqa E501
    assert o.note2 == 'T'
    t = astropy.time.Time(o.jdutc, format='jd').datetime
    assert t.year == 1988
    assert t.month == 10
    assert t.day == 5
    assert numpy.isclose(.046982, (t.hour + (t.minute + (t.second +
        t.microsecond / 1e6) / 60) / 60) / 24)
    h, m, s = astropy.coordinates.Angle(o.ra, unit='h').hms
    assert h == 0
    assert m == 52
    assert numpy.isclose(s, 27.557)
    sign, d, m, s = astropy.coordinates.Angle(o.dec, unit='deg').signed_dms
    assert sign == 1
    assert d == 38
    assert m == 16
    assert numpy.isclose(s, 5.91)
    assert o.mag == 10.74
    assert o.band == 'V'


def test_o10():
    """note2"""
    o = obs80.parseOpt('00433         C1993 09 29.84946 19 43 43.22 -14 20 51.0                z22620107')  # noqa E501
    assert o.note2 == 'C'


def test_satellite():
    o = obs80.parseSat('00433         S2016 05 18.90150 22 08 22.93 -13 53 27.3          16   RL~1tbGC51', '00433         s2016 05 18.90150 1 + 6328.2817 - 2371.9854 - 1216.2699   ~1tbGC51')  # noqa E501
    assert o.num == '00433'
    assert o.desig == ''
    assert o.disc == ''
    assert o.note1 == ''
    t = astropy.time.Time(o.jdutc, format='jd').datetime
    assert t.year == 2016
    assert t.month == 5
    assert t.day == 18
    assert numpy.isclose(.90150, (t.hour + (t.minute + (t.second +
        t.microsecond / 1e6) / 60) / 60) / 24)
    h, m, s = astropy.coordinates.Angle(o.ra, unit='h').hms
    assert h == 22
    assert m == 8
    assert numpy.isclose(s, 22.93)
    sign, d, m, s = astropy.coordinates.Angle(o.dec, unit='deg').signed_dms
    assert sign == -1
    assert d == 13
    assert m == 53
    assert numpy.isclose(s, 27.3)
    assert o.mag == 16
    assert o.band == 'R'
    assert o.cod == 'C51'
    assert numpy.isclose(o.x, 6328.2817 * 1e3 / astropy.constants.au.value)
    assert numpy.isclose(o.y, -2371.9854 * 1e3 / astropy.constants.au.value)
    assert numpy.isclose(o.z, -1216.2699 * 1e3 / astropy.constants.au.value)


def test_radar1():
    o = obs80.parseRad('00433         R1975 01 22.187500  15088536000                   430 251 JPLRS251', '00433         r1975 01 22.187500C        15000                      251 JPLRS251')  # noqa E501
    assert o.num == '00433'
    assert o.desig == ''
    assert o.note1 == ''
    t = astropy.time.Time(o.jdutc, format='jd').datetime
    assert t.year == 1975
    assert t.month == 1
    assert t.day == 22
    assert numpy.isclose(.187500, (t.hour + (t.minute + (t.second +
        t.microsecond / 1e6) / 60) / 60) / 24)
    assert o.delay == 150885360.00
    assert o.doppler is None
    assert o.freq == 430
    assert o.txcod == '251'
    assert o.rxcod == '251'
    assert o.ret == 'C'
    assert o.udelay == 15
    assert o.udoppler is None


def test_radar2():
    o = obs80.parseRad('00433         R1975 01 23.184028               -         130    430 251 JPLRS251', '00433         r1975 01 23.184028S                        2000       251 JPLRS251')  # noqa E501
    assert o.delay is None
    assert o.doppler == -1.3
    assert o.udelay is None
    assert o.udoppler == 2


def test_v():
    o = obs80.parseRov('     K03U11V  V2003 10 23.48631 01 37 13.60 +07 20 04.3          15.4 Rci8327247', '     K03U11V  v2003 10 23.48631 1 139.6119   +35.7803      70          ci8327247')  # noqa E501
    assert o.num == ''
    assert o.desig == 'K03U11V'
    assert o.disc == ''
    assert o.note1 == ''
    t = astropy.time.Time(o.jdutc, format='jd').datetime
    assert t.year == 2003
    assert t.month == 10
    assert t.day == 23
    assert numpy.isclose(.48631, (t.hour + (t.minute + (t.second +
        t.microsecond / 1e6) / 60) / 60) / 24)
    h, m, s = astropy.coordinates.Angle(o.ra, unit='h').hms
    assert h == 1
    assert m == 37
    assert numpy.isclose(s, 13.60)
    sign, d, m, s = astropy.coordinates.Angle(o.dec, unit='deg').signed_dms
    assert sign == 1
    assert d == 7
    assert m == 20
    assert numpy.isclose(s, 4.3)
    assert o.mag == 15.4
    assert o.band == 'R'
    assert o.lon == 139.6119
    assert o.lat == 35.7803
    assert o.alt == 70
