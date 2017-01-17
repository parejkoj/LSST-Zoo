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
set_lim=50
count=

#flag 1 is an asteroid in cutouts[i][3]
for i in range(len(cutouts)):
    if not np.isnan(np.sum(cutouts[i][0])) and cutouts[i][3]==0:
        count+=1
        """
        count+=1
        template=cutouts[i][0]
        science=cutouts[i][1]
        difference=cutouts[i][2]

        stitched_array=np.concatenate((template,science,difference),axis=1)
        minimum = stitched_array.min()
        Q=5
        scaled=luptonRGB.makeRGB(stitched_array,Q=Q,minimum=minimum)
        template=scaled[:,:20,0]
        science=scaled[:,20:40,0]
        difference=scaled[:,40:,0]
        cenx=9.5
        ceny=9.5

        f, (ax1, ax2, ax3) = plt.subplots(1, 3)

        ax1.imshow(template,cmap=cm.viridis,interpolation="none")
        ax1.axis('off')
        ax1.scatter(cenx, ceny, marker='+', color='r', s=100)
        ax2.imshow(science,cmap=cm.viridis,interpolation="none")
        ax2.axis('off')
        ax2.scatter(cenx, ceny, marker='+', color='r', s=100)
        ax3.imshow(difference,cmap=cm.viridis,interpolation="none")
        ax3.axis('off')
        ax3.scatter(cenx, ceny, marker='+', color='r', s=100)
        ax1.set_title("Template")
        ax2.set_title("Science")
        ax3.set_title("Difference")
        plt.tight_layout(w_pad=-3)
        plt.savefig('sub_sets/flagone/ast'+str(i)+'.jpg',bbox_inches='tight')
        plt.close()
     """
print count