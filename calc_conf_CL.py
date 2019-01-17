#!/usr/bin/env python

import numpy as np
import nibabel as nib
from dipy.tracking.metrics import downsample
from dipy.tracking.distances import bundles_distances_mdf
import matplotlib.pyplot as plt
from dipy.tracking.utils import length
import sys
import os


def makehist(x,name='none',saveloc=None, show=False):
    num_bins = 50
    # the histogram of the data
    fig, ax = plt.subplots(nrows=1, ncols=1)  # create figure & 1 axis
    #ax.plot([0,1,2], [10,20,3])
    n, bins, patches = ax.hist(x, num_bins, facecolor='green', alpha=0.5)
    # add a 'best fit' line
    #y = mlab.normpdf(bins, mu, sigma)
    #plt.plot(bins, y, 'r--')
    #plt.xlabel('Smarts')
    #plt.ylabel('Probability')
    ax.set_title(r'Histogram of Confidence: %s' % name)
    #fig.savefig(saveloc,dpi=100)
    #plt.close(fig)

    # Tweak spacing to prevent clipping of ylabel
    #plt.subplots_adjust(left=0.15)
    if show:
        plt.show()
    if saveloc:
        fig.savefig(saveloc)

def showsls(sls,values,outpath,show=False):
    from dipy.viz import window, actor,fvtk
    from dipy.data import fetch_bundles_2_subjects, read_bundles_2_subjects
    from dipy.tracking.streamline import transform_streamlines
    #renderer.clear()

    from dipy.tracking.streamline import length


    renderer=window.Renderer()

    hue = [0.5, 1]  # white to purple to red
    saturation = [0.0, 1.0]  # black to white

    lut_cmap = actor.colormap_lookup_table(
        scale_range=(values.min(), np.percentile(values,50)),
        hue_range=hue,
        saturation_range=saturation)

    stream_actor5 = actor.line(sls, values, linewidth=0.1,
                               lookup_colormap=lut_cmap)

    renderer.add(stream_actor5)
    bar3 = actor.scalar_bar(lut_cmap)

    renderer.add(bar3)

    # window.show(renderer, size=(600, 600), reset_camera=False)
    if outpath:
        window.record(renderer, out_path=outpath, size=(600, 600))
    if show:
        fvtk.show(renderer)

def doit(basedir, rootname, trkloc, debug=False):
    savename = '%s/%s_mtrx.txt' % (basedir, rootname)
    savetrk = '%s/%s.trk' % (basedir, rootname)
    savetrkss = '%s/%s_SS.trk' % (basedir,rootname)
    savetrkpic = savetrk.replace('.trk','.png')
    savetrkpicss = savetrk.replace('.trk','SS.png')
    savehist = '%s/%s_hist.png' % (basedir,rootname)
    subsampN = 12
    showhist=False
    showsls_var=False
    lenthr = 50
    pow = 1
    if debug:
        subsampN=12
        subset = 800

    trk,hdr = nib.trackvis.read(trkloc)
    sls = [item[0] for item in trk]
    lengths = list(length(sls))
    #print lengths
    print "cheerio"
    print len(lengths)
    print len(sls)
    print "two"
    print lengths[0]

    sls_long = []
    for n,j in enumerate(sls):
        if len(j) > 5:
            if lengths[n]>lenthr:
                print lengths[n]
                sls_long.append(j)
            else:
                print "BAD"
                print lengths[n]

    print "CHECKME"
    print len(sls)
    print len(sls_long)

    subsamp = [downsample(sl,subsampN) for sl in sls_long]
    if debug:
        subsamp = subsamp[0:subset]
        sls_long = sls_long[0:subset]
    score_mtrx = np.zeros([len(subsamp)])
    print "what"
    print len(subsamp)
    print len(sls)

    for i,sl in enumerate(subsamp):
        print str(i)+'/'+str(len(subsamp))
        mtrx2 = bundles_distances_mdf([sls_long[i]],sls_long)
        #mtrx2_oi = np.all([(mtrx2>0) & (mtrx2<5) & (not np.isnan(mtrx2))],axis=0)
        mtrx2_oi = (mtrx2>0) & (mtrx2<5) & ~np.isnan(mtrx2)
        zoom = mtrx2[mtrx2_oi]
        score = np.sum(np.divide(1,np.power(zoom,pow)))
        #makehist(zoom.T,score)
        score_mtrx[i] = score
        print score
        #print score_mtrx
    print savename
    print basedir
    print rootname
    np.savetxt(savename,score_mtrx)
    makehist(score_mtrx,name=rootname,saveloc=savehist,show=showhist)
    newhdr=hdr.copy()
    newhdr['n_properties']=1
    proplist=newhdr['property_name']
    proplist[0]='Cluster Confidence'
    newhdr['property_name']=proplist
    newtrkss = ((sl,None,score_mtrx[i]) for i,sl in enumerate(subsamp))
    newtrk = ((sl, None, score_mtrx[i]) for i,sl in enumerate(sls_long))
    nib.trackvis.write(savetrkss, newtrkss, newhdr)
    nib.trackvis.write(savetrk,newtrk,newhdr)
    #showsls(subsamp,score_mtrx,savetrkpicss,show=showsls_var)
    showsls(sls_long, score_mtrx, savetrkpic,show=showsls_var)


if __name__ == '__main__':
    trkloc = os.path.abspath(sys.argv[1])
    putloc = os.path.join(os.path.dirname(trkloc),'cci')
    if os.path.exists(trkloc):
        if not os.path.exists(putloc):
            os.mkdir(putloc)
        doit(putloc, os.path.basename(trkloc).replace('.trk','_cci'), trkloc, False)
