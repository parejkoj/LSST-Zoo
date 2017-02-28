import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
import os
import luptonRGB

import lsst.afw as afw
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.daf.persistence as dafPersist
from lsst.afw.math.warper import Warper

from esutil.sqlite_util import SqliteConnection

import lsst.log
logger = lsst.log.Log.getDefaultLogger()
logger.setLevel(lsst.log.WARN)

#
# Scary monkey patch!
#
from lsst.skymap.patchInfo import PatchInfo
def patch_hash(self):
    return hash(self._index)
PatchInfo.__hash__ = patch_hash


def zscale_image(input_img, contrast=0.25):
    """This emulates ds9's zscale feature. Returns the suggested minimum and
    maximum values to display."""

    samples = input_img.flatten()[::200]
    samples.sort()
    #if len(samples) < 60:
        #import pdb; pdb.set_trace()
    # Sort puts NaN values at the end, and
    # argmax returns the first NaN, so this avoids all NaNs.
    chop_start = int(0.10*np.argmax(samples))
    chop_end = int(0.90*np.argmax(samples))
    if chop_start == chop_end:
        return 1,1 #junk
    subset = samples[chop_start:chop_end]

    i_midpoint = int(len(subset)/2)
    #if len(subset) == 0:
        #import pdb; pdb.set_trace()
    I_mid = subset[i_midpoint]

    fit = np.polyfit(np.arange(len(subset)) - i_midpoint, subset, 1)
    # fit = [ slope, intercept]

    z1 = I_mid + fit[0]/contrast * (1-i_midpoint)/1.0
    z2 = I_mid + fit[0]/contrast * (len(subset)-i_midpoint)/1.0
    return z1, z2


def group_items(items, group_length):
    for n in xrange(0, len(items), group_length):
        yield items[n:(n+group_length)]


def make_patch_string(patch):
    return "{:d},{:d}".format(patch[0], patch[1])

def make_cutouts_for_patch(b, tract_info, patch_info, tracklets, cutout_size=30):

    visits = set()
    for s in tracklets['visits']:
        visits.update(s.split(","))

    visits = list(visits)

    wcs = tract_info.getWcs()
    dataref = {"tract": tract_info.getId(),
               "patch": make_patch_string(patch_info.getIndex()),
               "filter": "VR"}
    available_visits = filter(lambda v: butler.datasetExists("deepCoadd_tempExp",
                                                        visit=int(v),
                                                        **dataref),
                              visits)

    sci_cutouts =[]
    diff_cutouts=[]
    temp_cutouts=[]
    cutout_data=[]
    cached_matched_template = {}
    ccdnum = 10

    for visit in available_visits:
        for dataref in butler.subset('src', visit=int(visit), ccdnum=ccdnum):
            catalog = dataref.get('deepDiff_diaSrc')
            science=dataref.get('calexp', immediate=True) #Science Image
            difference=dataref.get('deepDiff_differenceExp', immediate=True) #difference image
            sci_wcs=science.getWcs()




            print("new science image")
            #import pdb; pdb.set_trace()
            #difference=dataref.get('deepDiff_differenceExp', immediate=True) #Difference
            count=0
            for source in catalog:
                SNR = source['base_PsfFlux_flux']/source['base_PsfFlux_fluxSigma']
                print SNR
                count+=1
                #import pdb; pdb.set_trace()
                scoord=source.getCoord() #[0] is RA [1] is DEC
                #import pdb; pdb.set_trace()
                #print scoord[0].asDegrees(), scoord[1].asDegrees()
                #sRA=source.get('coord_ra')
                #sDEC=source.get('coord_dec')

                #pixel_center = wcs.skyToPixel(afwCoord.IcrsCoord(scoord[0],scoord[1]))
                cutout_size=cutout_size
                #shifted_xy = pixel_center - source_patch.getOuterBBox().getBegin()
                #x_cut1=int((shifted_xy.getX() - cutout_size))
                #x_cut2=int((shifted_xy.getX() + cutout_size))
                #y_cut1=int((shifted_xy.getY() - cutout_size))
                #y_cut2=int((shifted_xy.getY() + cutout_size))

                sci_pixel_center=sci_wcs.skyToPixel(afwCoord.IcrsCoord(scoord[0],scoord[1]))
                scix_cut1=int((sci_pixel_center.getX() - cutout_size))
                scix_cut2=int((sci_pixel_center.getX() + cutout_size))
                sciy_cut1=int((sci_pixel_center.getY() - cutout_size))
                sciy_cut2=int((sci_pixel_center.getY() + cutout_size))

                if scix_cut1 < 0 or scix_cut2 < 0 or sciy_cut2 <0 or sciy_cut1 <0:
                    continue

                science_array=science.getMaskedImage().getImage().getArray()
                sci_cutout=science_array[sciy_cut1:sciy_cut2,scix_cut1:scix_cut2]
                #sci_cutout=np.fliplr(sci_cutout)
                #sci_cutout=np.flipud(sci_cutout)
                sci_cutouts.append(sci_cutout)

                difference_array=difference.getMaskedImage().getImage().getArray()
                diff_cutout=difference_array[sciy_cut1:sciy_cut2,scix_cut1:scix_cut2]
                #diff_cutout=np.fliplr(diff_cutout)
                #diff_cutout=np.flipud(diff_cutout)
                diff_cutouts.append(diff_cutout)


                source_patch=tract_info.findPatch(scoord)
                cache_key = (source_patch.getIndex(), int(visit), ccdnum)
                if cache_key in cached_matched_template.keys():
                    warpedTemplate = cached_matched_template[cache_key]
                else:
                    template=butler.get('deepCoadd_calexp',tract=tract_info.getId(),
                                        patch=make_patch_string(source_patch.getIndex()),
                                        filter="VR",immediate=True)

                    warper = Warper("lanczos4")
                    warpedTemplate = warper.warpExposure(science.getWcs(), template,
                                                         destBBox=science.getBBox())
                    cached_matched_template[cache_key] = warpedTemplate
                    print "regenerating " , cache_key

                template_image = warpedTemplate.getMaskedImage().getImage().getArray()
                temp_cutout=template_image[sciy_cut1:sciy_cut2,scix_cut1:scix_cut2]
                #import pdb; pdb.set_trace()
                temp_cutouts.append(temp_cutout)

                #Scaling: LuptonRGB
                #stitched_array=np.concatenate((temp_cutout,sci_cutout, diff_cutout),axis=1)
                #minimum = stitched_array.min()
                Q1=8
                Q2=15
                Q3=30
                #scaled=luptonRGB.makeRGB(stitched_array,Q=Q,minimum=minimum)
                #scaled_temp_cutout=scaled[:,:20,0]
                #scaled_sci_cutout=scaled[:,20:40,0]
                #scaled_diff_cutout=scaled[:,40:,0]
                #import pdb; pdb.set_trace()
                top_level_grid=gridspec.GridSpec(1,3)

                plt.subplot(top_level_grid[0]) #Template
                #scaled_temp_cutout=luptonRGB.makeRGB(temp_cutout,Q=Q,minimum=temp_cutout.min())
                z1,z2=zscale_image(sci_cutout)
                z1 /= 3.0
                z2 *= 3.0
                scaled_temp_cutout = ((temp_cutout - z1)/z2)
                plt.imshow(scaled_temp_cutout, interpolation="none", cmap=cm.viridis)
                plt.ylabel("Z-Scale", color="White")
                plt.title("Template",color="White")
                plt.axis("off")

                plt.subplot(top_level_grid[1]) #Science
                #scaled_sci_cutout = ((sci_cutout - z1)/z2) #.clip(0,1)
                #sci_min=sci_cutout.min()
                #scaled_sci_cutout=luptonRGB.makeRGB(sci_cutout,Q=Q,minimum=sci_cutout.min())
                z1,z2=zscale_image(sci_cutout)
                z1 /= 3.0
                z2 *= 3.0
                scaled_sci_cutout = ((sci_cutout - z1)/z2)
                plt.imshow(scaled_sci_cutout, interpolation="none", cmap=cm.viridis)
                plt.title("Science",color='white')
                plt.axis("off")

                plt.subplot(top_level_grid[2]) #Difference
                #scaled_sci_cutout = ((sci_cutout - z1)/z2) #.clip(0,1)
                #sci_min=sci_cutout.min()
                #scaled_diff_cutout=luptonRGB.makeRGB(diff_cutout,Q=Q,minimum=sci_cutout.min())
                z1,z2=zscale_image(diff_cutout)
                z1 /= 3.0
                z2 *= 3.0
                scaled_diff_cutout = ((diff_cutout - z1)/z2)
                plt.imshow(scaled_diff_cutout, interpolation="none", cmap=cm.viridis)
                plt.title("Difference",color='white')
                plt.axis("off")
                """
                #Lupton rows
                #Q1
                plt.subplot(top_level_grid[3])
                scaled_temp_cutout=luptonRGB.makeRGB(temp_cutout,Q=Q1,minimum=temp_cutout.min())
                plt.imshow(scaled_temp_cutout[:,:,0], interpolation="none", cmap=cm.viridis)
                plt.ylabel("Q = "+str(Q1), color="white")
                plt.axis("off")

                plt.subplot(top_level_grid[4])
                scaled_sci_cutout=luptonRGB.makeRGB(sci_cutout,Q=Q1,minimum=sci_cutout.min())
                plt.imshow(scaled_sci_cutout[:,:,0], interpolation="none", cmap=cm.viridis)
                plt.axis("off")

                plt.subplot(top_level_grid[5])
                scaled_diff_cutout=luptonRGB.makeRGB(diff_cutout,Q=Q1,minimum=sci_cutout.min())
                plt.imshow(scaled_diff_cutout[:,:,0], interpolation="none", cmap=cm.viridis)
                plt.axis("off")
                #Q2
                plt.subplot(top_level_grid[6])
                scaled_temp_cutout=luptonRGB.makeRGB(temp_cutout,Q=Q2,minimum=temp_cutout.min())
                plt.imshow(scaled_temp_cutout[:,:,0], interpolation="none", cmap=cm.viridis)
                plt.ylabel("Q = "+str(Q2))
                plt.axis("off")

                plt.subplot(top_level_grid[7])
                scaled_sci_cutout=luptonRGB.makeRGB(sci_cutout,Q=Q2,minimum=sci_cutout.min())
                plt.imshow(scaled_sci_cutout[:,:,0], interpolation="none", cmap=cm.viridis)
                plt.axis("off")

                plt.subplot(top_level_grid[8])
                scaled_diff_cutout=luptonRGB.makeRGB(diff_cutout,Q=Q2,minimum=sci_cutout.min())
                plt.imshow(scaled_diff_cutout[:,:,0], interpolation="none", cmap=cm.viridis)
                plt.axis("off")
                #Q3
                plt.subplot(top_level_grid[9])
                scaled_temp_cutout=luptonRGB.makeRGB(temp_cutout,Q=Q3,minimum=temp_cutout.min())
                plt.imshow(scaled_temp_cutout[:,:,0], interpolation="none", cmap=cm.viridis)
                plt.ylabel("Q = "+str(Q3), color="white")
                plt.axis("off")

                plt.subplot(top_level_grid[10])
                scaled_sci_cutout=luptonRGB.makeRGB(sci_cutout,Q=Q3,minimum=sci_cutout.min())
                plt.imshow(scaled_sci_cutout[:,:,0], interpolation="none", cmap=cm.viridis)
                plt.axis("off")

                plt.subplot(top_level_grid[11])
                scaled_diff_cutout=luptonRGB.makeRGB(diff_cutout,Q=Q3,minimum=sci_cutout.min())
                plt.imshow(scaled_diff_cutout[:,:,0], interpolation="none", cmap=cm.viridis)
                plt.axis("off")
                """


                plt.subplots_adjust(hspace=0.02, wspace=0.02, left=0.05, right=0.95,
                        top=0.95, bottom=0.05)
                plotfilename="cutouts/cutouts_tract{:d}_v{:d}_{:s}.png".format(tract_info.getId(),
                                                               int(visit),
                                                               str(count))
                print plotfilename
                plt.savefig(plotfilename,dpi=150, facecolor='k')

                #difference=difference.getMaskedImage().getImage().getArray()
                #diff_cutout=difference[x_cut1:x_cut2,y_cut1:y_cut2]
                #diff_cutouts.append(diff_cutout)
    return temp_cutouts, sci_cutouts, diff_cutouts

"""
#[5:48]
    patch_images = [butler.get("deepCoadd_tempExp", visit=int(v), **dataref) #Temporary exposure
                    for v in available_visits]

    patch_diffs = [butler.get("deepDiff_differenceExp", visit=int(v), ccd=1) #Temporary exposure
                    for v in available_visits]

    if len(patch_images) == 0:
        return [], [], (0,0)

    combined_img = np.sum(np.dstack([im.getMaskedImage().getImage().getArray() #Coadd
                                     for im in patch_images]), axis=-1)

    z1, z2 = zscale_image(combined_img)
    z1 /= 3.0
    z2 *= 3.0

    cutouts = []
    cutout_data = []
    sci_cutouts = []
    diff_cutouts=[]
    cuts=[]
    for tracklet in tracklets:

        # Some tracklets are generated from a different (but slightly
        # overlapping) tract, so they contain visits that aren't in
        # available_visits. Have to filter out these tracklets.
        this_tracklet_visits = tracklet['visits'].split(",")
        if  ((this_tracklet_visits[0] not in available_visits) and
             (this_tracklet_visits[1] not in available_visits)):
            continue
        #import pdb; pdb.set_trace()
        pixel_center = wcs.skyToPixel(afwCoord.IcrsCoord(tracklet['mean_ra']*afwGeom.degrees,
                                                         tracklet['mean_dec']*afwGeom.degrees))

        tracklet_length = (tracklet['delta_ra']*np.cos(np.degrees(tracklet['mean_dec'])) +
                           tracklet['delta_dec'])
        cutout_size = cutout_size

        shifted_xy = pixel_center - patch_info.getOuterBBox().getBegin()
        x_cut1=(shifted_xy.getX() - cutout_size)
        x_cut2=(shifted_xy.getX() + cutout_size)
        y_cut1=(shifted_xy.getY() - cutout_size)
        y_cut2=(shifted_xy.getY() + cutout_size)
        obj_cuts=[x_cut1,x_cut2,y_cut1,y_cut2]
        cuts.append(obj_cuts)
        cutout = combined_img[y_cut1:y_cut2,x_cut1:x_cut2]
        obj_sci_cutouts=[]
        for im in patch_images:
            print im
            im=im.getMaskedImage().getImage().getArray()
            sci_cutout=im[y_cut1:y_cut2,x_cut1:x_cut2]
            obj_sci_cutouts.append(sci_cutout)
        obj_diff_cutouts=[]
        for diff in patch_diffs:
            print diff
            diff=diff.getMaskedImage().getImage().getArray()
            diff_cutout=diff[y_cut1:y_cut2,x_cut1:x_cut2]
            obj_diff_cutouts.append(diff_cutout)
        sci_cutouts.append(obj_sci_cutouts)
        diff_cutouts.append(obj_diff_cutouts)
        cutouts.append(cutout)
        cutout_data.append({"tracklet_length": tracklet_length*3600.0})
    return cutouts, cutout_data, (z1, z2), sci_cutouts, diff_cutouts, cuts
"""

if __name__ == "__main__":
    path="/scratch/ctslater/NEO/"
    butler = dafPersist.Butler(os.path.join(path,"decam_NEO_repo/rerun/diffims1")) #butler instead of b
    db_conn = SqliteConnection("tracklets.db")
    db_conn.row_factory = None
    cursor = db_conn.cursor()

    skymap = butler.get("deepCoadd_skyMap")

    for target_tract in [0]:
        tract_info = skymap[target_tract]

        wcs = tract_info.getWcs()
        bbox = tract_info.getBBox()
        corner_coords = [wcs.pixelToSky(p + afwGeom.Extent2D(0,0)) for p in bbox.getCorners()]
        min_ra = min(c.getLongitude() for c in corner_coords).asDegrees()
        max_ra = max(c.getLongitude() for c in corner_coords).asDegrees()
        min_dec = min(c.getLatitude() for c in corner_coords).asDegrees()
        max_dec = max(c.getLatitude() for c in corner_coords).asDegrees()

        tracklets_res = cursor.execute("SELECT tracklet, avg(d.ra) as mean_ra, avg(d.dec) as mean_dec, "
                                       "max(d.ra) - min(d.ra) as delta_ra, "
                                       "max(d.dec) - min(d.dec) as delta_dec, "
                                       "group_concat(d.visit) as visits, count(*) "
                                       "FROM tracklet_links as tl JOIN detections as d "
                                       "ON tl.detection = d.id "
                                       "GROUP BY tracklet HAVING count(*) > 3 AND "
                                       "(mean_ra BETWEEN ? AND ?) AND (mean_dec BETWEEN ? AND ?)",
                                       (min_ra, max_ra, min_dec, max_dec))

        rows = tracklets_res.fetchall()
        dtype = db_conn._extract_row_dtype(cursor.description, rows)
        tracklets = np.array(rows, dtype=dtype)
        print("{:d} tracklets".format(len(tracklets)))

        tracklet_coords = [afw.coord.Coord(ra*afwGeom.degrees, dec*afwGeom.degrees)
                           for (ra, dec) in zip(tracklets['mean_ra'],
                                                tracklets['mean_dec'])]

        tracklet_coords_in_tract = filter(tract_info.contains, tracklet_coords)
        patches = [tract_info.findPatch(coord) for coord in tracklet_coords_in_tract]
        unique_patches = set(patches)

        for target_patch in list(unique_patches):
            if not butler.datasetExists("deepCoadd_calexp", tract=tract_info.getId(),
                          patch=make_patch_string(target_patch.getIndex()), filter="VR"):
                continue

            print("Patch ", target_patch)
            template=butler.get("deepCoadd_calexp", tract=tract_info.getId(), patch=make_patch_string(target_patch.getIndex()), filter="VR")
            template=template.getMaskedImage().getImage().getArray()
            #difference=butler.get("deepCoadd_calexp", tract=tract_info.getId(), patch=make_patch_string(target_patch.getIndex()), filter="VR")
            #difference=template.getMaskedImage().getImage().getArray()
            sel, = np.where(np.array(patches) == target_patch)
            cutouts, cutout_data, (z1, z2), sci_cutouts, diff_cutouts, cuts = make_cutouts_for_patch(butler, tract_info, target_patch, tracklets[sel])
            if len(cutouts) == 0:
                continue

            for group_n, cutout_group in enumerate(group_items(zip(cutouts, cutout_data, sci_cutouts, diff_cutouts, cuts), 1*1)): #4*4
                plt.figure(1).clear()
                #top_level_grid = gridspec.GridSpec(3, 2) #(4,4)
                print cuts
                for cutout_n, (cutout, cutout_data, sci_cutouts, diff_cutouts, cuts) in enumerate(cutout_group):
                    count=0
                    #print cutout_n
                    #print sci_cutouts
                    print cuts
                    template_cutout=template[cuts[2]:cuts[3],cuts[0]:cuts[1]]
                    try:
                        #temp_min=template_cutout.min()
                        #cut_min=cutout.min()
                        #Q=5
                        #scaled_template_cutout=luptonRGB.makeRGB(template_cutout,Q=Q,minimum=temp_min)
                        #scaled_cutout=luptonRGB.makeRGB(cutout,Q=Q,minimum=cut_min)
                        for sci_cutout, diff_cutout in sci_cutouts, diff_cutouts:
                            #print "made it to loop"
                            count+=1
                            top_level_grid=gridspec.GridSpec(1,4)

                            plt.subplot(top_level_grid[0]) #Template
                            scaled_template_cutout = ((template_cutout - z1)/z2).clip(0,1)
                            plt.imshow(scaled_template_cutout, interpolation="none", cmap=cm.viridis)
                            plt.title("Template",color='white')
                            plt.axis("off")

                            plt.subplot(top_level_grid[1]) #Science
                            scaled_sci_cutout = ((sci_cutout - z1)/z2) #.clip(0,1)
                            sci_min=sci_cutout.min()
                            #scaled_sci_cutout=luptonRGB.makeRGB(sci_cutout,Q=Q,minimum=sci_min)
                            plt.imshow(scaled_sci_cutout, interpolation="none", cmap=cm.viridis)
                            plt.title("Science",color='white')
                            plt.axis("off")

                            plt.subplot(top_level_grid[2]) #Difference
                            #diff_cutout=(template_cutout-sci_cutout)
                            #diff_min=difference.min()
                            #scaled_diff_cutout=luptonRGB.makeRGB(diff_cutout,Q=Q,minimum=diff_min)
                            scaled_diff_cutout = ((diff_cutout - z1)/z2).clip(0,1)
                            plt.imshow(scaled_diff_cutout, interpolation="none", cmap=cm.viridis)
                            plt.title("Difference",color='white')
                            plt.axis("off")

                            plt.subplot(top_level_grid[3]) #Coadd
                            scaled_cutout = ((cutout - z1)/z2).clip(0,1)
                            plt.imshow(scaled_cutout, interpolation="none", cmap=cm.viridis)
                            plt.title("Coadd", color='white')
                            plt.axis("off")

                            #print "made it to savefig"

                            plt.subplots_adjust(hspace=0.02, wspace=0.02, left=0.05, right=0.95,
                                    top=0.95, bottom=0.05)
                            plt.savefig("cutouts/cutouts_tract{:d}_p{:d}{:d}_{:0d}_{:s}.png".format(tract_info.getId(),
                                                                           target_patch.getIndex()[0],
                                                                           target_patch.getIndex()[1],
                                                                           group_n,str(count)),
                                                                           dpi=150, facecolor='k')
                            #print "past savefig"
                    except ValueError:
                        print "Bad Array"
