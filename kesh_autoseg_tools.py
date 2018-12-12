import nibabel as nib
import numpy as np

import re

from dipy.tracking.utils import length
from dipy.tracking import utils
from dipy.align.streamlinear import whole_brain_slr
from dipy.segment.clustering import qbx_and_merge
from dipy.segment.bundles import RecoBundles
from dipy.align.streamlinear import StreamlineLinearRegistration
from dipy.tracking.streamline import set_number_of_points

import nilearn.plotting as nip

from dipy.viz import window, actor


# Load nifti
def loadnii(niipath):
    im = nib.load(niipath)
    return im, im.get_data(), im.affine


# henry lab tracks saved with lowercase voxel order in header; convert it
def convert_trk_old2newapi(badfile):
    trk, hdr = nib.trackvis.read(badfile)
    new_hdr = hdr.copy()

    old_vox_order_UC = hdr['voxel_order'].item().decode('UTF-8').upper()
    new_hdr['voxel_order'] = old_vox_order_UC.encode('UTF-8')
    nib.trackvis.write(badfile.replace('.trk', '_newapi.trk'), trk, new_hdr)
    return badfile.replace('.trk', '_newapi.trk')


# load a tractogram with new nibabel streamline api
def loadtgm_newapi(trkpath):
    trkloaded = nib.streamlines.trk.TrkFile.load(trkpath)
    hdrloaded = trkloaded.header
    tg = trkloaded.tractogram
    return tg, hdrloaded


# filter by length
def filter_length(streamlines, minlen=40):
    print("calc lengths")
    lengths = list(length(streamlines))
    print("filter")
    long_sls = []
    for i, sl in enumerate(streamlines):
        if lengths[i] > minlen:
            long_sls.append(sl)
    return long_sls


# plot the roi
def plotroi(data, aff, bg=None):
    nip.plot_roi(nib.Nifti1Image(1*(data), aff), bg_img=bg, cmap='magma')


# filter based on freesurfer codes
def filter_freesurf(sls, apac_data, aff, code):
    return list(utils.target(sls, apac_data == code, affine=aff))


# make dictionary from excel sheet
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


# combine the rois into a single "NOT" exclusion and sets of "AND" inclusion
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


# generate a render
def genren_AGG(sls, sls2=None, niidata=None, roi1=None, roi2=None, roi3=None,
               aff=np.eye(4), putpath='test.png', showme=False,
               showaxes=False):

    renderer = window.Renderer()
    renderer.set_camera(position=(-606.93, -153.23, 28.70),
                        focal_point=(2.78, 11.06, 15.66),
                        view_up=(0, 0, 1))

    stream_actor = actor.line(sls)
    renderer.add(stream_actor)

    if sls2 is not None:
        stream_actor2 = actor.line(sls2, colors=(1, 1, 1))
        renderer.add(stream_actor2)

    if roi1 is not None:
        contour_actor1 = actor.contour_from_roi(roi1, affine=aff,
                                                color=(1., 1., 0.),
                                                opacity=0.5)
        renderer.add(contour_actor1)

    if roi2 is not None:
        contour_actor2 = actor.contour_from_roi(roi2, affine=aff,
                                                color=(1., 0., 0.),
                                                opacity=0.5)
        renderer.add(contour_actor2)
        
    if roi3 is not None:
        contour_actor3 = actor.contour_from_roi(roi3, affine=aff,
                                                color=(0., 0., 1.),
                                                opacity=0.5)
        renderer.add(contour_actor3)

    if niidata is not None:
        slice_actor = actor.slicer(niidata, affine=aff)
        renderer.add(slice_actor)

    if showaxes:
        axes = actor.axes()
        renderer.add(axes)

    if showme:
        window.show(renderer, size=(500, 500), reset_camera=False)
    window.record(renderer, out_path=putpath, size=(500, 500))
    # renderer.camera_info()
    del renderer
    return putpath


# targeting script to target streamlines with ROIs
def targetme(sls, include, exclude, aff):
    for i in range(include.shape[-1]):
        sls = list(utils.target(sls, include[:, :, :, i], affine=aff))
    for i in range(exclude.shape[-1]):
        sls = list(utils.target(sls, exclude[:, :, :, i], affine=aff,
                   include=False))
    return sls


# do a rough streamline registrion on the whole brain
def rough_reg(sub_fixed, temp_moving, save_base_abspath=None):
    # template moves to the subject space
    moved, transform, qb_cents1, qb_cents2 = whole_brain_slr(sub_fixed,
                                                             temp_moving,
                                                             verbose=True,
                                                             progressive=True)
    if save_base_abspath is not None:
        pklf = open(save_base_abspath+'.pkl', "wb")
        pickle.dump(transform, pklf)
    
    return moved, transform, qb_cents1, qb_cents2


# Recobundles wrapper
def run_rb(template, bucket, cluster_map=None, pruning_thr=10):
    # try pruning thresh 10 if not specific drop to 5
    if cluster_map is None:
        cluster_map = qbx_and_merge(bucket, thresholds=[40, 25, 20, 10])
    else:
        print("Loading provided cluster map")

    rb = RecoBundles(bucket, cluster_map=cluster_map, clust_thr=5)
    bundle_tsp, labels, bundle_bsp = rb.recognize(model_bundle=template,
                                                  model_clust_thr=5.,
                                                  reduction_thr=10,
                                                  pruning_thr=pruning_thr)
    return bundle_bsp, cluster_map


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


# streamline linear registration wrapper
def runslr(fixed, moving, npts=20, bounds=None, verbose=False):
    fixed_subsamp = set_number_of_points(fixed, npts)
    moving_subsamp = set_number_of_points(moving, npts)

    srr = StreamlineLinearRegistration(bounds=bounds, verbose=verbose)
    #print(srr.verbose)
    srm = srr.optimize(static=fixed_subsamp, moving=moving_subsamp)
    aligned = srm.transform(moving)
    return aligned


# replace this with your mapping between filenames and the excel column names
def convert_filename(mystring, file_list):
    mysplit = mystring.split(' ')
    hemi = mysplit[0][0]
    if mysplit[1] == 'SLF':
        if mysplit[-1][-1] == 'I':
            track = 'SLF.IP_'+mysplit[-1]
        elif mysplit[-1] == 'tp':
            track = 'SLF.tp'
    else:
        track = mysplit[-1]

    regex_string = '(?i)('+track+')_'+hemi+'.trk'
    prog = re.compile(regex_string)
    matches = []
    for i in file_list:
        if bool(prog.match(i)):
            matches.append(i)
    if len(matches) != 1:
        print('Need exactly one match for %s' % mystring)
        print(regex_string)
        print(matches)
        return None
    else:
        return matches[0]


# get the template and transform it into the space of the case with affine
def get_template(template_path, xfm_temp2case):
    template_tg, template_hdr = loadtgm_newapi(template_path)
    template_xfmd = template_tg.copy().apply_affine(xfm_temp2case).streamlines
    return template_xfmd
