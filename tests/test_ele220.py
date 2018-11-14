# mpcutilities/tests/test_ele220.py

"""
    
--------------------------------------------------------------

Nov 2018

Payne


Test the ele220 classes & functions written by Sonia Keys

These are for parsing & manipulating orbit data in MPC-220 format

Most / all of these tests were written by Sonia Keys 
 - Payne has done little-to-nothing to understand the 
 - breadth or depth of the testing 

--------------------------------------------------------------
    
"""

# Import third-party packages
# --------------------------------------------------------------
import numpy as np
import math

# Import the specific package/module/function we are testing
# --------------------------------------------------------------
import mpcutilities.ele220 as ele220


def test_to_julian_date():
    """
    Appropriately formatted description of what this test does ...
    """
    # test cases from Meeus
    assert ele220.to_julian_date(333, 1, 27.5) == 1842713.0
    assert ele220.to_julian_date(837, 4, 10.3) == 2026871.8
    assert ele220.to_julian_date(-123, 12, 31) == 1676496.5
    assert ele220.to_julian_date(-122, 1, 1) == 1676497.5
    assert ele220.to_julian_date(-1000, 7, 12.5) == 1356001
    assert ele220.to_julian_date(-1000, 2, 29) == 1355866.5
    assert ele220.to_julian_date(-1001, 8, 17.9) == 1355671.4
    assert ele220.to_julian_date(-4712, 1, 1.5) == 0
    # from mpc_astro.py
    assert ele220.to_julian_date(1957, 10, 4.81) == 2436116.31
    assert ele220.to_julian_date(2000, 1, 1.5) == 2451545.0
    assert ele220.to_julian_date(1582, 10, 4.9) == 2299160.4
    assert ele220.to_julian_date(1582, 10, 15.0) == 2299160.5


def test_unpackDate():
    # test cases from mpc_convert.py
    y, m, d, frac = ele220.unpackDate('K156R')
    assert y == 2015
    assert m == 6
    assert d == 27
    assert frac == 0.
    y, m, d, frac = ele220.unpackDate('-395P1179600')
    assert y == -239
    assert m == 5
    assert d == 25
    assert frac == .11796


def test_JDFromPacked():
    # test cases from mpc_convert.py
    j = ele220.JDFromPacked('K156R')
    assert j == 2457200.5
    j = ele220.JDFromPacked('-395P1179600')
    assert j == 1633907.61796


def test_220():
    e = ele220.Ele220('LM06C8L 0.15K16A4339018507176244313214879736009783459010378882420548564216      LM06C8L    19.7  K16BS  15.5  71.8 321.5   9.8 0.549   2.299   1     0     7 E XXXXXXXXXXMPC xxxxx K16C1 K16C1  0.22       0000             ')  # noqa E501
    assert e.desig() == 'LM06C8L'
    assert e.g() == .15
    assert e.timePeri() == 2457665.8390185
    assert e.argPeri() == 71.7624431
    assert e.node() == 321.4879736
    assert e.inc() == 9.783458999999999
    assert e.periDist() == 1.037888242
    assert e.ecc() == 0.548564216
    assert e.pertEpoch() is None
    assert e.pertScheme() == ''
    assert e.pertCoarse() == ''
    assert e.perturbers() == ''
    assert e.designation() == 'LM06C8L'
    assert e.h() == 19.7
    assert e.maEpoch() == 2457720.5
    assert e.ma1() == 15.5
    assert e.argPeri1() == 71.8
    assert e.node1() == 321.5
    assert e.inc1() == 9.8
    assert e.ecc3() == 0.549
    assert e.a3() == 2.299
    assert e.numOpp() == 1
    assert e.arc() == 0
    assert e.numObs() == 7
    assert e.u() == 10
    assert e.uNote() == 'E'
    assert e.comp() == 'XXXXXXXXX'
    assert e.ref() == 'MPC xxxxx'
    assert e.first() == 2457723.5
    assert e.last() == 2457723.5
    assert e.rms() == 0.22
    assert e.numScore() == -1
    assert e.curOppScore() == -1
    # compute full precision a
    k = 0.01720209895
    a = e.periDist() / (1 - e.ecc())
    assert abs(a - e.a3()) < 1e-3
    # compute full precision MA
    n = k / math.sqrt(a * a * a) * 180 / math.pi
    m = n * (e.maEpoch() - e.timePeri())
    if m < 0:
        m += 360
    assert abs(m - e.ma1()) < .1


def test_220_2():
    e = ele220.Ele220('K16P66L 0.15K165B812850902425636162311790513008660675311767228360146944491FK167V2016 PL66  21.5  K167V  48.8  24.3 231.2   8.7 0.147   1.379   1    10    15 6 MPCALB    MPO384740 K1686 K168G  0.27 h M-v 0038      4359236')  # noqa E501
    assert e.desig() == 'K16P66L'
    assert e.designation() == '2016 PL66'
    assert e.pertEpoch() == 2457600.5
    assert e.pertScheme() == 'DE403'
    assert e.pertCoarse() == 'M-v'
    assert e.perturbers() == '0038'
    assert e.numObs() == 15
    assert e.u() == 6
    assert e.uNote() == ''
    assert e.numScore() == 4
    assert e.curOppScore() == 3


def test_220_3():
    e = ele220.Ele220('00001   0.12K184R054656107281471130803142738010591701625585724940075705051 K167V     (1)    3.34 K167V 224.1  72.8  80.3  10.6 0.076   2.768 113 78700  6634 0 MPCLINUX  MPO384741 I011V K167N  0.60 h M-v 0030     z97    1')  # noqa E501
    assert e.desig() == '00001'
    assert e.designation() == '(1)'
    assert e.numObs() == 6634
    assert e.u() == 0
    assert e.uNote() == ''
    assert e.numScore() == 619
    assert e.curOppScore() == 7


def test_90():
    e = ele220.Ele90('K15C62F18.0 K14CJ394320818696268903071973369004844942218255187300177186443      K152FK152G')  # noqa E501
    assert e.desig() == 'K15C62F'
    assert e.h() == 18.0
    assert e.timePeri() == 2457010.8943208
    assert e.argPeri() == 186.96268899999998
    assert e.node() == 307.1973369
    assert e.inc() == 4.8449422
    assert e.periDist() == 1.8255187300000002
    assert e.ecc() == 0.177186443
    assert e.first() == 2457068.5
    assert e.last() == 2457069.5


def test_90_noH():
    e = ele220.Ele90('J26E00G     J263G826132008124575440770761267001043117824124847610112025420      J2635J263G')  # noqa E501
    assert e.h() is None


def test_sat():
    e = ele220.EleSat('    SK03J230K0544870259728524217390642831668148849592500922539820393086639 K047ES/2003 J 2316.7  K047E 275.0 285.2  64.3 148.8 0.393   0.152        29    16 0 MPCW      E2004-N18 K0326 K0337  0.35 h M-v 00385A    86    0')  # noqa E501
    assert e.permDesig() == ''
    assert e.tempDesig() == 'K03J230'
    assert e.pertEpoch() == 2453200.5
    assert e.designation() == 'S/2003 J 23'
    assert e.h() == 16.7
    assert e.maEpoch() == 2453200.5
    assert e.ma1() == 275.0
    assert e.argPeri1() == 285.2
    assert e.node1() == 64.3
    assert e.inc1() == 148.8
    assert e.ecc3() == 0.393
    assert e.a3() == 0.152
    assert e.arc() == 29
    assert e.numObs() == 16
    assert e.comp() == 'MPCW'
    assert e.ref() == 'E2004-N18'
    assert e.first() == 2452676.5
    assert e.last() == 2452705.5
    assert e.rms() == 0.35
    assert e.pertScheme() == 'DE403'
    assert e.pertCoarse() == 'M-v'
    assert e.perturbers() == '0038'
    assert e.numScore() == 8
    assert e.curOppScore() == 6
    e = ele220.EleSat('N010S       K143I381591014609319743195295901124377631202409545060223437421 K123EN X        11.0  K123E 132.2 146.1 319.5 124.4 0.223   0.310      2987    36 0 MPCW      MPC 78591 K018B K09AF  0.52 h M-v 00388A   255    0')  # noqa E501
    assert e.permDesig() == 'N010'
    assert e.tempDesig() == ''
    assert e.pertEpoch() == 2456000.5
    assert e.designation() == 'N X'
    assert e.h() == 11.0
    assert e.maEpoch() == 2456000.5
    assert e.ma1() == 132.2
    assert e.argPeri1() == 146.1
    assert e.node1() == 319.5
    assert e.inc1() == 124.4
    assert e.ecc3() == 0.223
    assert e.a3() == 0.310
    assert e.arc() == 2987
    assert e.numObs() == 36
    assert e.comp() == 'MPCW'
    assert e.ref() == 'MPC 78591'
    assert e.first() == 2452132.5
    assert e.last() == 2455119.5
    assert e.rms() == 0.52
    assert e.pertScheme() == 'DE403'
    assert e.pertCoarse() == 'M-v'
    assert e.perturbers() == '0038'
    assert e.numScore() == 25
    assert e.curOppScore() == 5


def test_cmt():
    e = ele220.EleComet('0001P-60K010-395P117960008811242880308093465163467390905853647000967587100a-3967Halley                       @YeoKia     1981MN    4.0 15   8.5 10  39 May  25.12    197,  643    161* G07   H59   h M-c 0008 9              2  +0.28    +0.0150               ')  # noqa E501
    assert e.permDesig() == '0001P'
    assert e.tempDesig() == 'P-60K010'
    assert e.frag() == ''
    assert e.pertEpoch() == 1633920.5
    assert e.pertScheme() == 'DE403'
    assert e.pertCoarse() == 'M-c'
    assert e.perturbers() == '0008'
    assert e.name() == 'Halley'
    assert e.comp() == '@YeoKia'
    assert e.pubYear() == '1981'
    assert e.pub() == 'MN'
    assert e.timeScale() == "UTC"
    assert e.pubVol() == '197,  643'
    assert e.numObs() == 161
    assert e.numObsAN() == ''
    assert e.numObsPlus() is False
    assert e.nonGrav() == ('2', 0.28, 0.015)
    assert e.first() == 2308004.5
    assert e.last() == 2363521.5
    assert e.notes() == ''
    assert e.rms() == -1


def test_cmtAN():
    e = ele220.EleComet('    CJ48L010J485F905235931705550362038209365023148926802076279080999874651aJ4859Honda-Bernasconi             Schrutka    1978QJ                     48 May  15.9052   19,   52  A 150  J4865 J4893       0000 8              0.000176 +0.000525 +0.001937 5    ')  # noqa E501
    assert e.permDesig() == ''
    assert e.tempDesig() == 'CJ48L010'
    assert e.frag() == ''
    assert e.pertEpoch() == 2432680.5
    assert e.pertScheme() == ''
    assert e.pertCoarse() == ''
    assert e.perturbers() == '0000'
    assert e.numObs() == 150
    assert e.numObsAN() == 'A'
    assert e.numObsPlus() is False


def test_cmtPlus():
    e = ele220.EleComet('    CJ47F020J4753919840018212855643230707517129155835709618359001000000000a     Becvar                       Cunningham  1947                       47 May   3.9198U HAC   807     3+  J4742 J474G       0000 0                                                ')  # noqa E501
    assert e.numObs() == 3
    assert e.numObsAN() == ''
    assert e.numObsPlus() is True
    assert e.notes() == ''
    assert e.rms() == -1
    assert e.pertEpoch() is None
    assert e.pertScheme() == ''
    assert e.pertCoarse() == ''
    assert e.perturbers() == ''


def test_cmtRMS():
    e = ele220.EleComet('0045P       K656R980508334482108360672994042012015530805916690530807694925 K657EHonda-Mrkos-Pajdusakova      MPCW        2011     13.5 20  19.5  5  11 Sept.28.7789  MPC 75294     84* K0144 K116C h M-v 003E 9  !         ! 2  +0.6644  -0.042993  +0.3656 0.6')  # noqa E501
    assert e.permDesig() == '0045P'
    assert e.tempDesig() == ''
    assert e.frag() == ''
    assert e.numObs() == 84
    assert e.notes() == '!'
    assert e.rms() == .6
    assert e.pertEpoch() == 2475480.5
    assert e.pertScheme() == 'DE403'
    assert e.pertCoarse() == 'M-v'
    assert e.perturbers() == '003E'


def test_cmtNotes():
    e = ele220.EleComet('    CI82R01bI829H724100606958529163476565632142011061200077506550999906817aI82A2Great September comet        Hufnagel    1919AN                     82 Sept.17.7241  209,   20  A1500  I8298 I835Q       0000 2 #R           0.000271 +0.012265 +0.012791 5    ')  # noqa E501
    assert e.frag() == 'b'
    assert e.notes() == "#R"
    assert e.pubYear() == '1919'
    assert e.pub() == 'AN'
    assert e.pubVol() == '209,   20'


def test_cmtFrag2():
    e = ele220.EleComet('0073P     aqK0667919538619876962900699107788011387323509392014510693302405      Schwassmann-Wachmann         MPCM        2006                       06 June  7.920   MPC 56799    228  K064E K0654       0000 0                                                ')  # noqa E501
    assert e.frag() == 'aq'


def test_cmtPub():
    e = ele220.EleComet('    CJ48V010J48AR427106910725206522110394996023116996401354209000999934980aJ48AGEclipse comet                Hirst       1954MNSA                   48 Oct. 27.4271U  13,   33    147  J48B8 J4943       0000 4              0.000032 +0.001294 +0.000518 6    ')  # noqa E501
    assert e.pubYear() == '1954'
    assert e.pub() == 'MNSA'
    assert e.pubVol() == '13,   33'


def test_cmtFirst():
    e = ele220.EleComet('    CE02D010E023L000000009100607971256950578054995792003800000001000000000a                                  Hind        1877NAT                    02 Mar. 21        16,   50         E021  E022        0000 0  Z                                             ')  # noqa E501
    assert e.first() == 2233138.5
    assert e.last() == 2233169.5
    assert e.nonGrav() is None


def test_cmtNonGrav():
    e = ele220.EleComet('0021P       K656F223909617141526591941230504031777208211002546680692613537 K6564Giacobini-Zinner             MPCW        2011      9.0 15  15.5  5  12 Feb. 11.7279  MPC 75010    917* J981K K1155 h M-c 0008 9  !         ! 2  +0.2542  -0.100893          0.8')  # noqa E501
    assert e.nonGrav() == ('2', .2542, -0.100893)
    e = ele220.EleComet('0045P       K656R980508334482108360672994042012015530805916690530807694925 K657EHonda-Mrkos-Pajdusakova      MPCW        2011     13.5 20  19.5  5  11 Sept.28.7789  MPC 75294     84* K0144 K116C h M-v 003E 9  !         ! 2  +0.6644  -0.042993  +0.3656 0.6')  # noqa E501
    assert e.nonGrav() == ('2', .6644, -.042993, .3656)


def test_cmtTimeScale():
    # old, default is UTC
    e = ele220.EleComet('    CI98U010I98BN652384012355204660977442659140344790007560109751000000000aI98B5Brooks                       Scharbe     1904AN                     98 Nov. 23.6524  164,  377    266  I98AM I98BQ       0000 3  E                                             ')  # noqa E501
    assert e.timeScale() == "UTC"
    # old, with TT/TDT specified
    e = ele220.EleComet('    CI99E010I994D478106000871005310264080958146269038703265758181000357024aI994GSwift                        Marsden     1978AJ                     99 Apr. 13.4781T  83,   64    124  I9935 I998B h M-c 0008 9              0.000009 -0.000109 -0.001253 7    ')  # noqa E501
    assert e.timeScale() == "TT/TDT"
    # new, default is TT/TDT
    e = ele220.EleComet('    CJ00B010J004T410827602440949670417942764146448611713315289221001058090aJ004BGiacobini                    Sekanina    1978AJ                     00 Apr. 29.4108   83,   64     32  J0023 J008I h M-c 0008 9              0.000030 +0.000057 -0.000722 7    ')  # noqa E501
    assert e.timeScale() == "TT/TDT"
    # new, with UTC specified
    e = ele220.EleComet('    CJ00O010J0083700477201242247143294136465062534237910148348011000410471aJ0089Borrelly-Brooks              @MS         1903EAN                    00 Aug.  3.7005U 1 (4)  16  A 400  J007O J00AR       0000 5              0.000105 +0.000610 -0.000565 5    ')  # noqa E501
    assert e.timeScale() == "UTC"


def test_cmt210():
    e = ele220.EleComet('0033P       K168M584299701907225590664897746022393861121602398720463027346 K1699Daniel                       MPCN        1997                       16 Aug. 22.5843  NK    656     38*    1979    -    1993    9  !          2  +0.3113  +0.074490          0.8')  # noqa E501
    assert e.first() is None
    assert e.last() is None
    assert e.pertEpoch() == 2457640.5
    assert e.pertScheme() == ''
    assert e.pertCoarse() == ''
    assert e.perturbers() == ''
    assert e.nonGrav() == ('2', 0.3113, 0.07449)


def test_fmt220():
    e0 = '10199   0.15K041J736304624244972373003998334023397900813.07603910172462412 K167V (10199)    6.6  K167V  71.8 242.4 300.4  23.4 0.172  15.801  18  9684   567 0 MPCLINUX  MPO365371 J88B5 K155C  0.45 h M-v 0038     B10    1'  # noqa E501
    e = ele220.Ele220(e0)
    o = ele220.EleFields(desig=e.desig(), g=e.g(), timePeri=e.timePeri(),
        argPeri=e.argPeri(), node=e.node(), inc=e.inc(), periDist=e.periDist(),
        ecc=e.ecc(), designation=e.designation(), h=e.h(), maEpoch=e.maEpoch(),
        numOpp=e.numOpp(), arc=e.arc(), numObs=e.numObs(), comp=e.comp(),
        ref=e.ref(), first=e.first(), last=e.last(), rms=e.rms())
    e1 = o.fmt220()
    print(e0)
    print(e1)
    assert len(e1) == 220
    assert e1[:e._ecc.stop] == e0[:e._ecc.stop]
    # undocumented field follows ecc.
    # then the source record has inexplicable justification for designation
    # so that is not tested here.
    assert e1[e._h.start:e._numObs.stop] == e0[e._h.start:e._numObs.stop]
    # then an undocumented field between numObs and comp, then the last
    # documented field is rms.
    assert e1[e._comp.start:e._rms.stop] == e0[e._comp.start:e._rms.stop]


def test_fmt90():
    e0 = 'J26A00C11.1 J261B801501008553655263325634812003351258918524740880139130865      J261AJ261B'
    e = ele220.Ele90(e0)
    o = ele220.EleFields(desig=e.desig(), h=e.h(), timePeri=e.timePeri(),
        argPeri=e.argPeri(), node=e.node(), inc=e.inc(), periDist=e.periDist(),
        ecc=e.ecc(), first=e.first(), last=e.last())
    assert o.fmt90() == e0
