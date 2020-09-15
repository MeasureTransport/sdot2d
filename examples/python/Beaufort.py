import pysdot as ot
import numpy as np
import matplotlib.pyplot as plt

import rasterio
from rasterio.features import shapes
import shapely.geometry as geom
#from shapely.geometry import shape, Polygon
from shapely.ops import unary_union
from descartes import PolygonPatch

from matplotlib.patches import Circle, Wedge, Patch
from matplotlib.collections import PatchCollection

from skimage import io
from skimage.color import rgb2gray


def ClipToIce(grayImg, lagDiag, xbnds, ybnds):

    polys = []
    for geo, value in shapes(threshImg, mask=(grayImg>0.7)):
        poly = geom.shape(geo).buffer(0) # Convert the rasterio result to a shapely shape
        polys.append(poly)

    icePoly = unary_union(polys).simplify(1.0,preserve_topology=True)

    Nx = float(grayImg.shape[1])
    Ny = float(grayImg.shape[0])

    clippedPolys = []
    for cellInd in range(lagDiag.NumCells()):
        verts = lagDiag.GetCellVertices(cellInd)

        pts = [( Nx*(verts[0,i]-xbnds[0])/(xbnds[1]-xbnds[0]), Ny*(verts[1,i]-ybnds[0])/(ybnds[1]-ybnds[0])) for i in range(verts.shape[1])]
        poly = geom.Polygon(pts)
        inter = icePoly.intersection(poly)
        if(type(inter)==geom.polygon.Polygon):
            clippedPolys.append(inter)
        elif(type(inter)==geom.multipolygon.MultiPolygon):
            clippedPolys += list(inter)

    # for each of the clipped polygons, make them convex
    convexPolys = []
    for poly in clippedPolys:
        pts = np.array(poly.exterior.coords.xy)
        otPoly = ot.Polygon(pts)
        parts = otPoly.ConvexPartition()
        convexPolys += [geom.Polygon(p.GetVertices().T) for p in ot.Polygon(pts).ConvexPartition()]

    return convexPolys, icePoly


thinBy = 1
img = io.imread('BeaufortLandfast_MODIS_04_04_2020.jpg')
img = img[600:900:thinBy, 1100:1400:thinBy]

threshold = 0.8
fig, axs = plt.subplots(ncols=2)
axs[0].imshow(img,cmap='Greys_r')

img = rgb2gray(img).astype('float32')
threshImg = ((img>threshold)+1e-3).astype("float32").T

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

opts = {'Lloyd Steps':200, 'Lloyd Tol':1e-3, 'Max Steps': 100, 'GTol Abs':1e-3}
diag = ot.SemidiscreteOT.BuildCentroidal(dist, numPts, opts)

clippedPolys, icePoly = ClipToIce(threshImg.T, diag, xbnds, ybnds)

# Plot the resulting centroidal Voronoi diagram
fig, axs = plt.subplots(ncols=3,figsize=(16,5))

axs[0].imshow(threshImg.T,extent=[xbnds[0],xbnds[1],ybnds[0],ybnds[1]],origin='lower',alpha=0.8,cmap='Greys_r')
axs[0].set_title('Thresholded Image')
axs[0].set_ylim(ybnds[1],ybnds[0])
axs[0].set_xlim(xbnds[0],xbnds[1])

ot.PlotDiagram(diag, axs[1], distribution=dist)
axs[1].set_title('Weighted CVD')
axs[1].set_ylim(ybnds[1],ybnds[0])


print('Number of clipped polynomials = ', len(clippedPolys))
patches = [PolygonPatch(p) for p in clippedPolys]
fig, ax = plt.subplots()
ax.add_collection(PatchCollection(patches,facecolor='gray',edgecolor='k',alpha=0.4))
ax.set_ylim(threshImg.shape[0],0)
ax.set_xlim(0,threshImg.shape[1])

plt.show()
