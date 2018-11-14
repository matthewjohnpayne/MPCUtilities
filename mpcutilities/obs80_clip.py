"""snip raw text observation records corresponding to selected fields,
operates on all plain text observation formats
Command line usage:
args: [datf] [{ field | - },...] 
datf -- raw observation text records (stdin)
field -- keyword in named obs80 tuple objects selecting values to be output
  "-" means pass through the input data record in this position; see -a
opts: 
#-a -- all fields (in place of giving all the fields as command arguments)

ISB 8/2017
"""

import sys, os
# package parent directory must be in sys.path
from collections import namedtuple
from components.obs80.obs80 import *
from components.obs80 import obs80


def parseOpt(line):
    n2 = obs80.note2(line)
    if n2 not in 'PeCTMcEOHNn AXx':
        raise ValueError('note2 = "{}"'.format(n2))

    return line, Optical(
        num=line[:5].strip(),
        desig=line[5:12].strip(),
        disc=line[12].strip(),
        note1=line[13].strip(),
        note2=n2.strip(),
        jdutc=line[15:32].strip(),
        ra=line[32:44].strip(),
        dec=line[44:56].strip(),
        mag=line[65:70].strip(),
        band=line[70].strip(),
        cod=line[77:80]
    )

obs80.parseOpt = parseOpt

au_km = obs80.au_km
floatPx = obs80.floatPx

def parseSat(line1, line2):
    """ returns Spacebased for successful parse, raises exception otherwise"""
    
    ax = floatPx(line2[34:45].strip())
    ay = floatPx(line2[46:57].strip())
    az = floatPx(line2[58:69].strip())
    if line2[32] == '1':
        ax /= au_km
        ay /= au_km
        az /= au_km
    
    origs = (line1+line2).replace('\n',r'\n')
    return origs, SpaceBased(
        # line1
        num=line1[:5].strip(),
        desig=line1[5:12].strip(),
        disc=line1[12].strip(),
        note1=line1[13].strip(),
        jdutc=line1[15:32].strip(),
        ra=line1[32:44].strip(),
        dec=line1[44:56].strip(),
        mag=line1[65:70].strip(),
        band=line1[70].strip(),
        cod=line1[77:80].strip(),
        # line2
        x=str(ax),
        y=str(ay),
        z=str(az)
    )

obs80.parseSat = parseSat

def parseRoving(line1, line2):
    """ returns Roving for successful parse, raises exception otherwise"""
    
    origs = (line1+line2).replace('\n',r'\n')
    return origs, Roving(
        # line1
        num=line1[:5].strip(),
        desig=line1[5:12].strip(),
        disc=line1[12].strip(),
        note1=line1[13].strip(),
        jdutc=line1[15:32].strip(),
        ra=line1[32:44].strip(),
        dec=line1[44:56].strip(),
        mag=line1[65:70].strip(),
        band=line1[70].strip(),
        # line2
        lon=line2[34:44].strip(),
        lat=str(floatPx(line2[45:55])),
        alt=line2[56:61].strip()
    )

obs80.parseRov = parseRoving

def parseRadar(line1, line2):
    """ returns Radar for successful parse, raises exception otherwise"""
    
    origs = (line1+line2).replace('\n',r'\n')
    return origs, Radar(
        # line1
        num=line1[:5].strip(),
        desig=line1[5:12].strip(),
        note1=line1[13].strip(),
        jdutc=line1[15:32].strip(),
        delay=line1[32:43] + '.' + line1[43:47] if line1[42] != ' ' else None,
        doppler=str(floatPx(
            line1[47:58]+'.'+line1[58:62])) if line1[57] != ' ' else None,
        freq=line1[62:67] + '.' + line1[67:68] + line2[62:68],
        txcod=line1[68:71],
        rxcod=line1[77:80],
        # line2
        ret=line2[32],
        udelay=line2[33:43]+'.'+line2[43:47] if line2[42] != ' ' else None,
        udoppler=line2[47:58]+'.'+line2[58:62] if line2[57] != ' ' else None
    )

obs80.parseRad = parseRadar

def split_obs(strm):
    for s in strm:
        eol = s.find(r'\n')
        if eol > 0:
            s = s.rstrip().replace(r'\n','\n').replace('\\\\','\\')
            yield s[:eol+1]
            if eol+1 < len(s):
                yield s[eol+1:]
        else:
            yield s

def run(src, flds):
    """print CSV formatted input for Postgres COPY FROM"""
    def qq(s): 
        s = s.replace('"','""').replace("'","''")
        return s if s.find(',') < 0 else '"%s"' % s
    for t in parse80(split_obs(src)):
        if not isinstance(t[0],Exception):
            s, obs = t
            print(','.join([getattr(obs,k,"") if k != '-' else qq(
                s.rstrip()) for k in flds]))
        else:
            print('parse80: Error, {}'.format(t), file=sys.stderr)


def main():
    src = len(sys.argv) > 1 and sys.argv[1]
    if not (src and os.path.exists(src)):
        src = sys.stdin
    else:
        src = open(src)
        sys.argv.pop(0)
    flds = sys.argv[1].split(',') if len(sys.argv) > 1 else [
        'cod','jdutc','ra','dec']
    run(src, flds)


if __name__ == '__main__':
    main()
