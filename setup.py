from distutils.core import setup, Extension

import os, glob, numpy, sys

def indir(dir, files): return [dir+f for f in files]

setup(name='pfits',
    ext_modules = [
        Extension('pfits',
            ['src/pfits_wrap.cpp'] + \
            indir('src/cfitsio/', [
'buffers.c', 'cfileio.c', 'checksum.c', 'compress.c', 'drvrfile.c', 'drvrmem.c', 'drvrnet.c', 'drvrsmem.c', 'drvrgsiftp.c', 'editcol.c', 'edithdu.c', 'eval_l.c', 'eval_y.c', 'eval_f.c', 'fitscore.c', 'getcol.c', 'getcolb.c', 'getcold.c', 'getcole.c', 'getcoli.c', 'getcolj.c', 'getcolk.c', 'getcoll.c', 'getcols.c', 'getcolsb.c', 'getcoluk.c', 'getcolui.c', 'getcoluj.c', 'getkey.c', 'group.c', 'grparser.c', 'histo.c', 'iraffits.c', 'modkey.c', 'putcol.c', 'putcolb.c', 'putcold.c', 'putcole.c', 'putcoli.c', 'putcolj.c', 'putcolk.c', 'putcoluk.c', 'putcoll.c', 'putcols.c', 'putcolsb.c', 'putcolu.c', 'putcolui.c', 'putcoluj.c', 'putkey.c', 'region.c', 'scalnull.c', 'swapproc.c', 'wcssub.c', 'wcsutil.c', 'imcompress.c', 'quantize.c', 'ricecomp.c', 'pliocomp.c', 'fits_hcompress.c', 'fits_hdecompress.c',]),
            include_dirs = [numpy.get_include(), 'src/cfitsio'],
        )
    ],
)
            
