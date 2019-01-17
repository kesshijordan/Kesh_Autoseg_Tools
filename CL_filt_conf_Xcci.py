#!/usr/bin/env python

import numpy as np
import nibabel as nib
from dipy.tracking.metrics import downsample
from dipy.tracking.distances import bundles_distances_mdf
import matplotlib.pyplot as plt
from dipy.tracking.utils import length
import sys
import os

#in the quickbundles paper, the authors removed streamlines less than 40mm,
#which seems like as good a threshold as any


def doit(trkloc, ccithr, temphdr):
    if ccithr != int(ccithr):
        ccithr_str = 'pt'.join(str(ccithr).split('.'))
    else:
        ccithr_str = str(ccithr)
    suffix = '_ccithr'+ccithr_str+'.trk'
    savetrk = trkloc.replace('.trk', suffix)

    trk,hdr = nib.trackvis.read(trkloc)
    sls = [item[0] for item in trk]
    cci = [item[2][0] for item in trk]

    sls_good = []
    for n,j in enumerate(sls):
        if cci[n] > ccithr:
            sls_good.append(j)

    print("CHECKME")
    print(len(sls))
    print(len(sls_good))

    newtrk = ((sl, None, None) for sl in sls_good)
    nib.trackvis.write(savetrk,newtrk,temphdr)


if __name__ == '__main__':
    trkloc = os.path.abspath(sys.argv[1])
    ccithr = float(sys.argv[2])
    template = sys.argv[3]

    ttrk,temphdr = nib.trackvis.read(template)
    doit(trkloc, ccithr, temphdr)
