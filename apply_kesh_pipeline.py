import pandas as pd
import os
import dipy as dp
from dipy.tracking.utils import length
from dipy.tracking import utils
from glob import glob
import nibabel as nib
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

import nilearn.plotting as nip

from dipy.viz import window, actor, ui
from IPython.display import Image

from dipy.align.streamlinear import whole_brain_slr, slr_with_qb, transform_streamlines
from dipy.segment.clustering import qbx_and_merge
from dipy.segment.bundles import RecoBundles

track_of_interest = 'Left ARCUATE'

template_basepath = './qb_templates/Base_CTRL'
wb_template = os.path.join(
    template_basepath, 'Whole_Brain_long_resaved_newapi.trk')

template_boi = '/Users/kesshijordan/Desktop/IU_Bloomington/qb_templates/Arcuate_template.trk'

putdir = './pyAFQuicobundles_freesurfer_test_outputs'
if not os.path.exists(putdir):
    os.mkdir(putdir)

roi_matrix_path = 'WM_pathways_freesurfer_based_seg_EC_KMJmod.xlsx'

pathpath = './pyafquicobundles_freesurfer_test_basepath.txt'


# Load files
apac_path = glob(os.path.join(case_basepath, 'aparc+aseg_di*.nii*'))[0]
fa_path = glob(os.path.join(case_basepath, '*_fa.nii*'))[0]
md_path = glob(os.path.join(case_basepath, '*_md.nii*'))[0]
wb_path = glob(os.path.join(case_basepath, 'Whole_Brain*newapi.trk'))[0]

with open(pathpath, 'r') as myfile:
    case_basepath = myfile.read().replace('\n', '').replace('\'', '')

wb_template_tg, wb_template_hdr = loadtgm_newapi(wb_template)
wb_template_sls = wb_template_tg.streamlines
template_arcuate_tg, template_arcuate_hdr = loadtgm_newapi(template_arcuate)
template_arcuate_sls = template_arcuate_tg.streamlines


def loadnii(niipath):
    im = nib.load(niipath)
    return im, im.get_data(), im.affine


def loadtgm_newapi(trkpath):
    trkloaded = nib.streamlines.trk.TrkFile.load(trkpath)
    hdrloaded = trkloaded.header
    tg = trkloaded.tractogram
    return tg, hdrloaded


def filter_length(streamlines, minlen=40):
    print("calc lengths")
    lengths = list(length(streamlines))
    print("filter")
    long_sls = []
    for i, sl in enumerate(streamlines):
        if lengths[i] > minlen:
            long_sls.append(sl)
    return long_sls


def filter_freesurf(sls, apac_data, aff, code):
    return list(utils.target(sls, apac_data == code, affine=aff))


apac_im, apac_data, apac_aff = loadnii(apac_path)
fa_im, fa_data, fa_aff = loadnii(fa_path)
md_im, md_data, md_aff = loadnii(md_path)

roi_matrix = pd.read_excel(roi_matrix_path)


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


mydict = build_dict(track_of_interest, roi_matrix=roi_matrix)
include, exclude = combine_rois(mydict, apac_data)


def targetme(sls, include, exclude, aff):
    for i in range(include.shape[-1]):
        sls = list(utils.target(sls, include[:, :, :, i], affine=aff))
    for i in range(exclude.shape[-1]):
        sls = list(utils.target(
            sls, exclude[:, :, :, i], affine=aff, include=False))
    return sls


bundle = targetme(long_sls, include, exclude, fa_aff)


def rough_reg(sub_fixed, temp_moving):
    # template moves to the subject space
    # qb_thr=5 errored
    moved, transform, qb_centroids1, qb_centroids2 = whole_brain_slr(sub_fixed, temp_moving,
                                                                     verbose=True,
                                                                     progressive=True)
    return moved, transform, qb_centroids1, qb_centroids2


def run_rb(templatesls, bucketosls, cluster_map=None, pruning_thr=10):
    # try pruning thresh 10 if not specific drop to 5
    if cluster_map is None:
        cluster_map = qbx_and_merge(bucketosls, thresholds=[40, 25, 20, 10])
    else:
        print("Loading provided cluster map")

    rb = RecoBundles(bucketosls, cluster_map=cluster_map, clust_thr=5)
    recognized_atlassp, rec_labels, recognized_ptsp = rb.recognize(model_bundle=templatesls,
                                                                   model_clust_thr=5.,
                                                                   reduction_thr=10, pruning_thr=pruning_thr)
    '''rb = RecoBundles(bucketosls, cluster_map=cluster_map, clust_thr=10)
    recognized, rec_labels, rec_trans = rb.recognize(model_bundle=templatesls,
                                                         model_clust_thr=1.)'''
    #D = bundles_distances_mam(templatesls, recognized)

    return recognized_ptsp, cluster_map


from dipy.align.streamlinear import StreamlineLinearRegistration
from dipy.tracking.streamline import set_number_of_points


def runslr(fixed, moving, npts=20):
    fixed_subsamp = set_number_of_points(fixed, npts)
    moving_subsamp = set_number_of_points(moving, npts)

    srr = StreamlineLinearRegistration()
    srm = srr.optimize(static=fixed_subsamp, moving=moving_subsamp)
    aligned = srm.transform(moving)
    return aligned


moved_temp2case, xfm_temp2case, qbc1_temp2pcase, qbc2_temp2case = rough_reg(
    long_sls[::10], wb_template_sls)
template_boi_sls_xfmd = template_boi_tg.copy(
).apply_affine(xfm_temp2case).streamlines

template_boi_sls_xfmd_slr = runslr(arc, template_boi_sls_xfmd, npts=10)

template_arrayseq = nib.streamlines.array_sequence.ArraySequence(
    template_boi_sls_xfmd)
bucket_arrayseq = nib.streamlines.array_sequence.ArraySequence(boi)

rb_boi, ig_ip = run_rb(template_arrayseq, bucket_arrayseq)
