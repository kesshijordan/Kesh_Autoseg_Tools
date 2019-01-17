#!/usr/bin/env python

'''
This script takes two .trk files and combines them into one track file. 
python CL_merge_trks.py track1_path.trk track2_path.trk merged_path.trk 
'''

import numpy as np
import nibabel as nib
import sys
import os


def doit(trk1path,trk2path,putpath):
    trk,hdr = nib.trackvis.read(trk1path)
    sls = [item[0] for item in trk]
    
    trk2,hdr2 = nib.trackvis.read(trk2path)
    sls2 = [item[0] for item in trk2]

    newhdr = hdr.copy()
    newhdr['n_properties']=1
    proplist=newhdr['property_name']
    proplist[0] = 'id'
    newhdr['property_name']=proplist

    newtrk = list(((sl, None, np.array(0)) for sl in sls))
    newtrk2 = list(((sl, None, np.array(1)) for sl in sls2))
    nib.trackvis.write(putpath,newtrk+newtrk2,newhdr)
    print(len(newtrk), len(newtrk2), len(newtrk+newtrk2))

if __name__ == '__main__':
    trkloc1 = os.path.abspath(sys.argv[1])
    trkloc2 = os.path.abspath(sys.argv[2])
    putpath = os.path.abspath(sys.argv[3])

    doit(trkloc1, trkloc2, putpath)
