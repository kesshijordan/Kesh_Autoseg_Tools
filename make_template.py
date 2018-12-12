from glob import glob
from IPython.display import Image

import numpy as np
from nibabel import trackvis as tv
from dipy.segment.clustering import QuickBundles
from dipy.data import get_data
from dipy.viz import fvtk
import matplotlib as mpl
import seaborn as sns
import pickle

from dipy.tracking.streamline import length
from dipy.align.streamlinear import StreamlineLinearRegistration
#from streamlinear_kmjmod import StreamlineLinearRegistration
from dipy.tracking.streamline import set_number_of_points
import os

import sys

import nibabel as nib

from dipy.viz import window, actor
from dipy.viz.fvtk import camera

from dipy.tracking.distances import bundles_distances_mdf
from dipy.tracking.utils import move_streamlines

from dipy.align.streamlinear import whole_brain_slr

import kesh_autoseg_tools as kat

# Path to a whole-brain streamline dataset in MNI space
# https://figshare.com/articles/Simple_model_bundle_atlas_for_RecoBundles/6483614
MNI_wholebrain_trkfile_path = '/Users/kesshijordan/ref_data/Atlas_in_MNI_Space_16_bundles/whole_brain/whole_brain_MNI.trk'

# The Keystone subject will be used to map between your cohort and MNI space
Keystone_subject_wholebrain_trkfile_path = '/Users/kesshijordan/Desktop/IU_Bloomington/qb_templates/Base_CTRL/Whole_Brain_long_resaved_newapi.trk'
# Also include the bundle of interest (BOI) segmented in the Keystone subject
Keystone_subject_BOI_trkfile_path = '/Users/kesshijordan/Desktop/IU_Bloomington/qb_templates/Base_CTRL/Results/Arcuate_L.trk'

# These are the bundles that will be included in the template
BOI_list = glob(
    '/Users/kesshijordan/ref_data/control_manual_data/*/Arcuate_L.trk')
BOI_list

mni_tg, mni_hdr = kat.loadtgm_newapi(MNI_wholebrain_trkfile_path)
mni_sls = mni_tg.streamlines

key_tg, key_hdr = kat.loadtgm_newapi(Keystone_subject_wholebrain_trkfile_path)
key_sls = key_tg.streamlines

keystone_inMNI, transform_key2mni, qb_cents1, qb_cents2 = kat.rough_reg(
    mni_sls, key_sls)

key_boi_tg, key_boi_hdr = kat.loadtgm_newapi(Keystone_subject_BOI_trkfile_path)
key_boi_sls = key_boi_tg.streamlines


def load_trk_list(path_list):
    bundle_list = []
    for p in path_list:
        tg, hdr = kat.loadtgm_newapi(p)
        sls = tg.streamlines
        bundle_list.append(sls)
    return bundle_list


bundle_list = load_trk_list(BOI_list)


def make_kesh_template(bundle_list, keystone_boi, qb_thresh=5., Nsubsamp=20,
                       clsz_thresh=5, keystone2MNI_xfm=None, verbose=False):
    '''
    bundle_list: list of independent bundles (lists) not assumed to be in the same space
    keystone_boi: bundle (list) of streamlines that will be the anchor bundle all
                  others are registered to for the template
    qb_thresh: threshold for quickbundle (determines how finely each bundle is clustered)
    Nsubsamp: subsampling for quickbundles and SLR
    clsz_thresh: how many streamlines a cluster must have to be included in the template*
    keystone2MNI_SLR: streamlinear registration between the whole brain keystone and MNI**
    verbose: if you want to print info about each bundle as it runs set this to True


    *qb_thresh adn clsz_thresh are related. If you have a fine parcellation
    (low qb_thresh) then the clsz_threshold should be quite low since clusters
    will be small.

    **PROVIDE THIS IF (and only if) YOU WANT THE RESULT TO BE IN MNI SPACE OTHERWISE
    IT WILL BE IN KEYSTONE SPACE
    '''

    kesh_template_sls = []
    rejected_sls = []
    boi_sls_subsamp = set_number_of_points(keystone_boi, Nsubsamp)
    for i, sls in enumerate(bundle_list):
        print(len(bundle_list)-i)
        sls_subsamp = set_number_of_points(sls, Nsubsamp)
        qb = QuickBundles(threshold=qb_thresh)
        clusters = qb.cluster(sls)
        cluster_sizes = [len(cl) for cl in clusters]
        # enforce that clusters smaller than a threshold are not in template
        centroids = clusters.centroids
        slr = StreamlineLinearRegistration()
        srm = slr.optimize(static=boi_sls_subsamp, moving=sls_subsamp)
        xfmd_centroids = srm.transform(centroids)
        # NOTE: we actually want to upsample the centroids so the template has
        # better properties... what's the most efficient way to do that?
        for j, b in enumerate(xfmd_centroids):
            if cluster_sizes[j] < clsz_thresh:
                rejected_sls.append(xfmd_centroids.pop(j))
        kesh_template_sls += xfmd_centroids
        if verbose:
            print('Bundle %i' % i)
            print('N centroids: %i' % len(centroids))
            print('kept %i rejected %i total %i' % (len(kesh_template_sls),
                                                    len(rejected_sls), len(clusters)))
    if keystone2MNI_xfm:
        print('MNI YAY!')
    return kesh_template_sls, rejected_sls


template_sls, rejected_sls = make_kesh_template(bundle_list, key_boi_sls)

template_sls_mni = list(move_streamlines(template_sls, transform_key2mni))

template_tgm_mni = nib.streamlines.tractogram.Tractogram(
    streamlines=template_sls_mni, affine_to_rasmm=mni_tg.affine_to_rasmm)

nib.streamlines.save(template_tgm_mni, 'testing_arcuate_template_raw.trk')
