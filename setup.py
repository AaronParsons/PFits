from distutils.core import setup, Extension

import os, glob, numpy, sys
if 'upload' in sys.argv or 'register' in sys.argv:
    from ez_setup import use_setuptools; use_setuptools()
    from setuptools import setup, Extension

__version__ = open('VERSION').read().strip()

def get_description():
    lines = [L.strip() for L in open('README').readlines()]
    d_start = None
    for cnt, L in enumerate(lines):
        if L.startswith('DESCRIPTION'): d_start = cnt + 1
        elif not d_start is None:
            if len(L) == 0: return ' '.join(lines[d_start:cnt])
    raise RuntimeError('Bad README')

def indir(path, files):
    return [os.path.join(path, f) for f in files]

setup(name='pfits',
    version = __version__,
    description = 'A Python FITS interface built using CFITSIO',
    long_description = get_description(),
    license = 'GPL',
    author = 'Aaron Parsons',
    author_email = 'aparsons@astron.berkeley.edu',
    url = 'http://pypi.python.org/pypi/pfits',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    setup_requires = ['numpy>=1.2'],
    ext_modules = [
        Extension('pfits',
            ['src/pfits_wrap.cpp'] + \
            indir('src/cfitsio/', [
'buffers.c', 'cfileio.c', 'checksum.c', 'compress.c', 'drvrfile.c', 'drvrmem.c', 'drvrnet.c', 'drvrsmem.c', 'drvrgsiftp.c', 'editcol.c', 'edithdu.c', 'eval_l.c', 'eval_y.c', 'eval_f.c', 'fitscore.c', 'getcol.c', 'getcolb.c', 'getcold.c', 'getcole.c', 'getcoli.c', 'getcolj.c', 'getcolk.c', 'getcoll.c', 'getcols.c', 'getcolsb.c', 'getcoluk.c', 'getcolui.c', 'getcoluj.c', 'getkey.c', 'group.c', 'grparser.c', 'histo.c', 'iraffits.c', 'modkey.c', 'putcol.c', 'putcolb.c', 'putcold.c', 'putcole.c', 'putcoli.c', 'putcolj.c', 'putcolk.c', 'putcoluk.c', 'putcoll.c', 'putcols.c', 'putcolsb.c', 'putcolu.c', 'putcolui.c', 'putcoluj.c', 'putkey.c', 'region.c', 'scalnull.c', 'swapproc.c', 'wcssub.c', 'wcsutil.c', 'imcompress.c', 'quantize.c', 'ricecomp.c', 'pliocomp.c', 'fits_hcompress.c', 'fits_hdecompress.c',]),
            include_dirs = [numpy.get_include(), 'src/cfitsio'],
        )
    ],
)
            
