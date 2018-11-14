import io, sys
from . import obs80

namedtuple = obs80.collections.namedtuple
Roving = namedtuple('Roving', list(obs80.Roving._fields) + ['cod'])


def is_ob(x):
    """filter for valid 80-char observation records"""
    return not isinstance(x[0], Exception)


def gen_obs(s):
    for o in filter(is_ob, obs80.parse80(s)):
        if o.__class__ == obs80.Roving:
            kw = dict([(f, getattr(o, f)) for f in o._fields])
            kw['cod'] = '247'
            o = Roving(**kw)
        yield o


def scan_obs80(text):
    """wrapper for obs80.parse80 that cleans and buffers the input
    observation text, returns generator of named tuples"""
    def clean(s):
        s = s.replace("''", "'").replace("\\\\", "\\")
        s = s.rstrip().replace(r"\n", "\n").replace("\n\n", "\n")
        return s if s and s[-1] == '\n' else s + '\n'
    B = io.StringIO()
    for s in text:
        B.write(clean(s))
    B.seek(0)
    return gen_obs(B)
