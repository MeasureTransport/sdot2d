import pysdot as ot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import rasterio
from rasterio.features import shapes
import shapely.geometry as geom
from shapely.ops import unary_union
from descartes import PolygonPatch
from skimage import io
from skimage.color import rgb2gray
from scipy.stats import powerlaw
import time

def ClipToIce(grayImg, lagDiag, xbnds, ybnds):

    # Make ice extent polygon for clipping diagram cells
    polys = []
    for geo, value in shapes(grayImg, mask=(grayImg>0.7)):
        poly = geom.shape(geo).buffer(0) # Convert the rasterio result to a shapely shape
        polys.append(poly)

    icePoly = unary_union(polys).simplify(1.0,preserve_topology=True)

    Nx = float(grayImg.shape[1])
    Ny = float(grayImg.shape[0])

    allPolys = []

    # Clip polygons to ice extent and save boundary polygons
    clippedPolys = []

    # Polygons entirely within ice extent - not clipped by boundaries
    insidePolys = []

    # Polygons clipped by boundaries
    edgePolys = []

    for cellInd in range(lagDiag.NumCells()):
        verts = lagDiag.GetCellVertices(cellInd)

        pts = [( Nx*(verts[0,i]-xbnds[0])/(xbnds[1]-xbnds[0]), Ny*(verts[1,i]-ybnds[0])/(ybnds[1]-ybnds[0])) for i in range(verts.shape[1])]
        poly = geom.Polygon(pts)

        allPolys.append(poly)

        # Check for particles inside ice extent, but not crossing the boundary
        if(icePoly.intersects(poly)):
            if(icePoly.overlaps(poly)==False):
                insidePolys.append(poly)
            else:
                edgePolys.append(poly)

        # Save polygons clipped by ice extent
        inter = icePoly.intersection(poly)
        if(type(inter)==geom.polygon.Polygon):
            clippedPolys.append(inter)
        elif(type(inter)==geom.multipolygon.MultiPolygon):
            clippedPolys += list(inter)

    # Make sure each clipped polygon is convex
    convexPolys = []
    for poly in clippedPolys:
        pts = np.array(poly.exterior.coords[:])
        if (pts.size > 0):
            # convexPolys += [geom.Polygon(p.GetVertices().T) for p in ot.Polygon(pts.T).ConvexPartition()]
            p = ot.Polygon(pts.T).ConvexHull()
            convexPolys += [geom.Polygon(p.GetVertices().T)]

    return convexPolys, clippedPolys, insidePolys, edgePolys, allPolys

def GenerateSeeds(numPts, img, grid, tol=1e-2):
    """ Generates seed points that lie on regions of thresholded image with
        intensities greater than the specified tolerance `tol`.  This is useful
        for seeding points over the ice before constructing a centroidal
        Laguerre diagram.

        Uses simple rejection sampling to randomly place points inside the ice region.
    """

    pts = np.zeros((2,numPts))
    for i in range(numPts):

        # Just keep drawing random samples until we find one over the ice
        numAttempts = 0
        maxAttempts = 100000
        while(numAttempts<maxAttempts):
            numAttempts += 1
            # Generate a point in the bounding box
            samp = np.array([ grid.xMin + (grid.xMax-grid.xMin)*np.random.rand(1),
                              grid.yMin + (grid.yMax-grid.yMin)*np.random.rand(1)])

            # Get the index of the pixel containing this point
            xInd = grid.LeftNode(samp[0])
            yInd = grid.BottomNode(samp[1])

            # Check if the thresholded image is nonzero at this location
            if(img[xInd,yInd]>tol):
                pts[:,i] = samp[:,0]
                break

        assert(numAttempts<maxAttempts)

    return pts

# Start timer
tic = time.perf_counter()

# Read in image, thin and threshold
thinSchedule = [40,20,1]
threshold = 0.5

# Number of polygons desired in final diagram
numPts = 1000
print("Number of seed points: {0}".format(numPts))


full_img = io.imread('NaresStrait_06_26_2003_Trimmed.jpg')
full_img = rgb2gray(full_img).astype('float32').T

# Create grid and bounding box
xbnds = [0.0,full_img.shape[0]/np.min(full_img.shape)] # minimum and maximum x values
ybnds = [0.0,full_img.shape[1]/np.min(full_img.shape)] # minimum and maximum y values
bbox = ot.BoundingBox(xbnds[0],xbnds[1],ybnds[0],ybnds[1])

full_grid = ot.RegularGrid(bbox, full_img.shape[0], full_img.shape[1])

# Random generate intial seed points
seedPts = GenerateSeeds(numPts,full_img,full_grid,threshold)
seedProbs = np.ones(numPts) / numPts

# Use centroidal diagrams with coarse images to inform finer solves
for thinBy in thinSchedule:

    print('=========================================')
    print('Working on thin level {}'.format(thinBy))
    print('-----------------------------------------\n')

    img = full_img[::thinBy, ::thinBy]
    threshImg = ((img>threshold)).astype("float32")

    # Construct the centroidal diagram with capacity constraints
    grid = ot.RegularGrid(bbox, threshImg.shape[0], threshImg.shape[1])

    threshImg /= np.sum(threshImg*grid.dx*grid.dy)
    dist = ot.DiscretizedDistribution(grid, threshImg)

    # Construct the Centroidal Voronoi diagram
    opts = {'Lloyd Steps':200, 'Lloyd Tol':1e-8, 'Max Steps': 100, 'GTol Abs':1e-4*np.sqrt(numPts)}
    diag = ot.SemidiscreteOT.BuildCentroidal(dist, seedPts, seedProbs, opts)

    # Update the seed points for the next iteration
    seedPts = diag.Centroids(dist)

    # Grab areas of each cell
    areas = diag.Areas(dist)

    # Plot the resulting centroidal Voronoi diagram
    fig, axs = plt.subplots(ncols=2,figsize=(16,5))

    axs[0].imshow(threshImg.T,extent=[xbnds[0],xbnds[1],ybnds[0],ybnds[1]],origin='lower',alpha=0.8,cmap='Greys_r')
    axs[0].set_title('Thresholded Image (thinBy={})'.format(thinBy))
    axs[0].set_ylim(ybnds[1],ybnds[0])
    axs[0].set_xlim(xbnds[0],xbnds[1])

    ot.PlotDiagram(diag, axs[1], distribution=dist, cell_colors=areas)
    axs[1].set_title('Weighted CVD (thinBy={})'.format(thinBy))
    axs[1].set_ylim(ybnds[1],ybnds[0])

    plt.savefig('NaresStrait_Centroidal_{:03d}.png'.format(thinBy))
    plt.close()



# Domain bounds
convexPolys, clippedPolys, insidePolys, edgePolys, allPolys = ClipToIce(threshImg.T, diag, xbnds, ybnds)

# End timer
toc = time.perf_counter()
print("\nProcess took {0:.4f} seconds".format(toc - tic))

# Plot the resulting centroidal Voronoi diagram
fig, axs = plt.subplots(ncols=2,figsize=(16,5))

axs[0].imshow(threshImg.T,extent=[xbnds[0],xbnds[1],ybnds[0],ybnds[1]],origin='lower',alpha=0.8,cmap='Greys_r')
axs[0].set_title('Thresholded Image')
axs[0].set_ylim(ybnds[1],ybnds[0])
axs[0].set_xlim(xbnds[0],xbnds[1])

ot.PlotDiagram(diag, axs[1], distribution=dist, cell_colors=areas)
axs[1].set_title('Weighted CVD')
axs[1].set_ylim(ybnds[1],ybnds[0])

plt.show()
