# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 12:31:27 2016

@author: Doug
"""

import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from matplotlib import cm
import luptonRGB

data=scipy.io.loadmat("cutouts_threedet_tracklets_303482")
cutouts=data["cutouts"]
Q=5

for i in range(len(cutouts)):
    if not np.isnan(np.sum(cutouts[i][0])):
        template=cutouts[i][0]
        science=cutouts[i][1]
        difference=cutouts[i][2]

        stitched_array=np.concatenate((template,science,difference),axis=1)
        minimum = stitched_array.min()
        Q=5
        scaled=luptonRGB.makeRGB(stitched_array,Q=Q,minimum=minimum)
        template=scaled[:,:20,0]
        science=scaled[:,21:41,0]
        difference=scaled[:,40:,0]

        f, (ax1, ax2, ax3) = plt.subplots(1, 3)

        ax1.imshow(template,cmap=cm.viridis,interpolation="none")
        ax1.axis('off')
        ax2.imshow(science,cmap=cm.viridis,interpolation="none")
        ax2.axis('off')
        ax3.imshow(difference,cmap=cm.viridis,interpolation="none")
        ax3.axis('off')
        ax1.set_title("Template")
        ax2.set_title("Science")
        ax3.set_title("Difference")
        plt.tight_layout(w_pad=-1.5)
        plt.savefig('sub_sets/Asteroids/ast'+str(i)+'.jpg',bbox_inches='tight')
        plt.close()