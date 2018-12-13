#!/usr/bin/env python
import nibabel as nib
import sys
'''
This is a quick script to address an incompatability between how .trk 
files have been saved historically (in the henry lab) and the new nibabel 
API. The loading does not accept lowercase voxel order, and the choice was 
not to extend the API to lowercase (instead enforce the file format), which 
means historical .trk files need to be converted.
'''
badfile = sys.argv[1]

trk,hdr = nib.trackvis.read(badfile)
new_hdr = hdr.copy()
new_hdr['voxel_order'] = new_hdr['voxel_order'].item().decode('UTF-8').upper().encode('UTF-8')
nib.trackvis.write(badfile.replace('.trk','_newapi.trk'),trk,new_hdr)

