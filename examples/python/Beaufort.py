import pysdot as ot
import numpy as np
import matplotlib.pyplot as plt

from skimage import io
from skimage.color import rgb2gray

img = io.imread('BeaufortLandfast_MODIS_04_04_2020.jpg')

thinBy = 1
threshold = 0.7
fig, axs = plt.subplots(ncols=2)
axs[0].imshow(img[600:900:thinBy, 1100:1400:thinBy],cmap='Greys_r')

img = rgb2gray(img)
threshImg = ((img[600:900:thinBy, 1100:1400:thinBy]>threshold)+1e-3).astype("float").T
axs[1].imshow(threshImg.T,cmap='Greys_r')
axs[0].set_title('Original Image')
axs[1].set_title('Thresholded Image')


numPts = 1000
xbnds = [0.0,threshImg.shape[0]/np.min(threshImg.shape)] # minimum and maximum x values
ybnds = [0.0,threshImg.shape[1]/np.min(threshImg.shape)] # minimum and maximum y values

bbox = ot.BoundingBox(xbnds[0],xbnds[1],ybnds[0],ybnds[1])
grid = ot.RegularGrid(bbox, threshImg.shape[0], threshImg.shape[1])

# Construct the centroidal diagram with capacity constraints
dist = ot.DiscretizedDistribution(grid,threshImg)

diag = ot.SemidiscreteOT.BuildCentroidal(dist, numPts,50, 0.003)

# Plot the resulting centroidal Voronoi diagram
fig, axs = plt.subplots(ncols=2,figsize=(14,6))
ot.PlotDiagram(diag, axs[0], distribution=dist)
axs[0].set_title('Weighted CVD')
axs[0].set_ylim(ybnds[1],ybnds[0])

axs[1].imshow(threshImg.T,extent=[xbnds[0],xbnds[1],ybnds[0],ybnds[1]],origin='lower',alpha=0.8,cmap='Greys_r')
axs[1].set_ylim(ybnds[1],ybnds[0])
axs[1].set_title('Thresholded Image')

plt.show()
