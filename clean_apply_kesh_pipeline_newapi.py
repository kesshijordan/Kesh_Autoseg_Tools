#!/usr/bin/python

'''
This runs a merged freesurfer targeting/recobundles approach that both constrains endpoints and matches shape. Run this by providing a patient wholebrain streamline dataset .trk file, a patient aparc+aseg.nii from freesurfer, and an output directory. If you're in the ALBA or Henry Labs at UCSF, then make sure
'''

import pandas as pd
import os

from dipy.tracking import utils

import nibabel as nib
import numpy as np


from dipy.align.streamlinear import whole_brain_slr
from dipy.segment.bundles import RecoBundles
from dipy.io.streamline import load_trk


def loadnii(niipath):
    im = nib.load(niipath)
    return im, im.get_data(), im.affine


# load a tractogram with new nibabel streamline api
def loadtgm_newapi(trkpath):
    trkloaded = nib.streamlines.trk.TrkFile.load(trkpath)
    hdrloaded = trkloaded.header
    tg = trkloaded.tractogram
    return tg, hdrloaded


# save with the old track format
def save_old_trk(streamlines, hdr, img, savepath):

    aff_rasmm2tv = nib.streamlines.trk.get_affine_rasmm_to_trackvis(hdr)

    new_hdr = hdr.copy()
    new_hdr['voxel_order'] = "".join(nib.orientations.aff2axcodes(img.affine))
    new_hdr['dimensions'] = img.get_data().shape
    new_hdr['voxel_sizes'] = img.header.get_zooms()
    new_hdr['voxel_to_rasmm'] = img.affine.copy()

    save_tgm = nib.streamlines.tractogram.Tractogram(streamlines=streamlines,
                                                     affine_to_rasmm=np.eye(4))

    save_tgm_xfmd = save_tgm.copy().apply_affine(aff_rasmm2tv)

    save_trk = nib.streamlines.TrkFile(save_tgm_xfmd, header=new_hdr)

    nib.streamlines.save(save_trk, filename=savepath)


# henry lab tracks saved with lowercase voxel order in header; convert it
def convert_trk_old2newapi(badfile):
    trk, hdr = nib.trackvis.read(badfile)
    new_hdr = hdr.copy()

    old_vox_order_UC = hdr['voxel_order'].item().decode('UTF-8').upper()
    new_hdr['voxel_order'] = old_vox_order_UC.encode('UTF-8')
    nib.trackvis.write(badfile.replace('.trk', '_newapi.trk'), trk, new_hdr)
    return badfile.replace('.trk', '_newapi.trk')


def filter_freesurf(sls, apac_data, aff, code):
    return list(utils.target(sls, apac_data == code, affine=aff))


def build_dict(track_name, roi_matrix):
    rois = roi_matrix[track_name]
    roi_groups = list(rois.value_counts().index)

    roi_dict = {}
    roi_dict['include'] = {}
    roi_dict['exclude'] = {}

    for i in roi_groups:
        if i > 0:
            roitype = 'include'
        elif i < 0:
            roitype = 'exclude'
        else:
            raise('UHOH')
        setname = 'set'+str(int(i))
        roi_dict[roitype][setname] = {}

        temp = roi_matrix[rois == i]
        for j, name in enumerate(temp['VOIS']):
            roi_dict[roitype][setname][name] = temp['aparc+aseg'].values[j]
    return roi_dict


def combine_rois(mydict, apac):
    x, y, z = apac.shape
    include = np.zeros([x, y, z, len(mydict['include'])])
    exclude = np.zeros([x, y, z, len(mydict['exclude'])])
    for i, iset in enumerate(mydict['include'].keys()):
        for item in mydict['include'][iset].items():
            include[:, :, :, i] += 1*(apac == item[-1])
    for j, jset in enumerate(mydict['exclude'].keys()):
        for item in mydict['exclude'][jset].items():
            exclude[:, :, :, j] += 1*(apac == item[-1])
    return include, exclude


def targetme(sls, include, exclude, aff):
    for i in range(include.shape[-1]):
        sls = list(utils.target(sls, include[:, :, :, i], affine=aff))
    for i in range(exclude.shape[-1]):
        sls = list(utils.target(
            sls, exclude[:, :, :, i], affine=aff, include=False))
    return sls


def rough_reg(sub_fixed, temp_moving):
    # template moves to the subject space
    # qb_thr=5 errored
    moved, transform, qb_centroids1, qb_centroids2 = whole_brain_slr(
        sub_fixed, temp_moving, verbose=True, progressive=True)
    return moved, transform, qb_centroids1, qb_centroids2


# Recobundles wrapper
def run_rb(template, bucket, pruning_thr=10):
    # try pruning thresh 10 if not specific drop to 5

    rb = RecoBundles(bucket, clust_thr=5)
    # TODO: for efficiency, we want to segment all model bundes at once
    recog_bundle, recog_labels = rb.recognize(model_bundle=template,
                                              model_clust_thr=5.,
                                              reduction_thr=10,
                                              pruning_thr=pruning_thr,
                                              reduction_distance='mam')
    return recog_bundle, recog_labels


def runkeshpype(wb_path, apac_path, putdir, convertapi=False):

    # This is the whole brain streamline dataset in the template space
    wb_template = '/Users/kesshijordan/Desktop/cleanup_desktop2/IU_Bloomington/qb_templates/Base_CTRL/Whole_Brain_long_resaved_newapi.trk'

    # This is the lookup table for what freesurfer ROIs are targeted per track
    roi_matrix_path = 'WM_pathways_freesurfer_based_seg_EC_KMJmod.xlsx'

    # This is the lookup table for what bundle templates are used per track
    readme_path = 'kesh_pipeline_bundle_template_lookup.xlsx'
    readme_df = pd.read_excel(readme_path).dropna(subset=['file_name'])

    # Build a dictionary of template paths
    template_path_dict = {}
    for c in readme_df.column_name:
        template_path_dict[c] = readme_df[readme_df.column_name ==
                                          c].file_name.values[0]
    print('Segmenting the following:')
    print(template_path_dict)

    if not os.path.exists(putdir):
        print('Making output directory')
        print(putdir)
        os.mkdir(putdir)

    if convertapi:
        wb_path = convert_trk_old2newapi(wb_path)

    print('loading whole brain of case')
    wb_tg, wb_hdr = loadtgm_newapi(wb_path)
    wb_sls = wb_tg.streamlines

    print('loading whole brain of template')
    wb_template_sls, wb_template_hdr = load_trk(wb_template)

    apac_im, apac_data, apac_aff = loadnii(apac_path)

    roi_matrix = pd.read_excel(roi_matrix_path)

    moved_temp2case, xfm_temp2case, qbc1_temp2pcase, qbc2_temp2case = rough_reg(
        wb_sls, wb_template_sls)

    # BUNDLE SPECIFIC CODE
    for track_of_interest in template_path_dict.keys():
        print('Segmenting %s' % track_of_interest)
        mydict = build_dict(track_of_interest, roi_matrix=roi_matrix)
        include, exclude = combine_rois(mydict, apac_data)

        template_boi = template_path_dict[track_of_interest]
        template_boi_sls, template_boi_hdr = load_trk(template_boi)
        template_boi_tg = nib.streamlines.Tractogram(template_boi_sls)
        bundle = targetme(wb_sls, include, exclude, apac_aff)
        template_boi_sls_xfmd = template_boi_tg.copy(
        ).apply_affine(xfm_temp2case).streamlines

        template_arrayseq = nib.streamlines.array_sequence.ArraySequence(
            template_boi_sls_xfmd)
        bucket_arrayseq = nib.streamlines.array_sequence.ArraySequence(bundle)
        rb_boi, ig_ip = run_rb(template_arrayseq, bucket_arrayseq)

        savepath = os.path.join(
            putdir, track_of_interest.replace(' ', '_')+'_TESTING.trk')
        save_old_trk(bucket_arrayseq[ig_ip], wb_hdr,
                     apac_im, savepath)
        print('Saved to %s' % savepath)


if __name__ == '__main__':
    import sys
    case_whole_brain = sys.argv[1]
    case_aparc_aseg = sys.argv[2]
    case_putpath = sys.argv[3]

    runkeshpype(case_whole_brain, case_aparc_aseg, case_putpath,
                convertapi=True)
