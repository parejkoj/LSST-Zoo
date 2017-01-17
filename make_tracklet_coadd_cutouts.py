import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
import os

import lsst.afw as afw
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.daf.persistence as dafPersist

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
    if len(samples) < 60:
        import pdb; pdb.set_trace()
    # Sort puts NaN values at the end, and
    # argmax returns the first NaN, so this avoids all NaNs.
    chop_start = int(0.10*np.argmax(samples))
    chop_end = int(0.90*np.argmax(samples))
    if chop_start == chop_end:
        return 1,1 #junk
    subset = samples[chop_start:chop_end]

    i_midpoint = int(len(subset)/2)
    if len(subset) == 0:
        import pdb; pdb.set_trace()
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
    patch_images = [butler.get("deepCoadd_tempExp", visit=int(v), **dataref) #Coadd?
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
    for tracklet in tracklets:

        # Some tracklets are generated from a different (but slightly
        # overlapping) tract, so they contain visits that aren't in
        # available_visits. Have to filter out these tracklets.
        this_tracklet_visits = tracklet['visits'].split(",")
        if  ((this_tracklet_visits[0] not in available_visits) and
             (this_tracklet_visits[1] not in available_visits)):
            continue

        pixel_center = wcs.skyToPixel(afwCoord.IcrsCoord(tracklet['mean_ra']*afwGeom.degrees,
                                                         tracklet['mean_dec']*afwGeom.degrees))

        tracklet_length = (tracklet['delta_ra']*np.cos(np.degrees(tracklet['mean_dec'])) +
                           tracklet['delta_dec'])
        cutout_size = cutout_size

        shifted_xy = pixel_center - patch_info.getOuterBBox().getBegin()
        cutout = combined_img[(shifted_xy.getY() - cutout_size):(shifted_xy.getY() + cutout_size),
                              (shifted_xy.getX() - cutout_size):(shifted_xy.getX() + cutout_size)]
        cutouts.append(cutout)
        cutout_data.append({"tracklet_length": tracklet_length*3600.0})
    return cutouts, cutout_data, (z1, z2)


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

            sel, = np.where(np.array(patches) == target_patch)
            cutouts, cutout_data, (z1, z2) = make_cutouts_for_patch(butler, tract_info, target_patch, tracklets[sel])
            if len(cutouts) == 0:
                continue

            for group_n, cutout_group in enumerate(group_items(zip(cutouts, cutout_data), 1*1)): #4*4
                plt.figure(1).clear()
                top_level_grid = gridspec.GridSpec(1, 1) #(4,4)

                for cutout_n, (cutout, cutout_data) in enumerate(cutout_group):
                    plt.subplot(top_level_grid[cutout_n])
                    scaled_cutout = ((cutout - z1)/z2).clip(0,1)
                    plt.imshow(scaled_cutout, interpolation="none" cmap=cm.viridis)
                    plt.title("cutouts_tract{:d}_p{:d}{:d}_{:0d}.png".format(tract_info.getId(),
                                                                               target_patch.getIndex()[0],
                                                                               target_patch.getIndex()[1],
                                                                               group_n"))
                    plt.axis("off")

                plt.subplots_adjust(hspace=0.02, wspace=0.02, left=0.05, right=0.95,
                                    top=0.95, bottom=0.05)
                plt.savefig("cutouts/cutouts_tract{:d}_p{:d}{:d}_{:0d}.png".format(tract_info.getId(),
                                                                           target_patch.getIndex()[0],
                                                                           target_patch.getIndex()[1],
                                                                           group_n),
                            dpi=150, facecolor='k')
