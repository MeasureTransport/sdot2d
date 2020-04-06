import pysdot as ot
import numpy as np
import matplotlib.pyplot as plt



numPts = 100

xbnds = [0.0,1.0] # minimum and maximum x values
ybnds = [0.0,1.0] # minimum and maximum y values

bbox = ot.BoundingBox(xbnds[0],xbnds[1],ybnds[0],ybnds[1])

# Construct the Centroidal Voronoi diagram.  This function uses Lloyd's algorithm
# with latin hypercube samples as initial points (https://en.wikipedia.org/wiki/Lloyd%27s_algorithm)
# Arguments to BuildCentroidal are:
# - The bounding box
# - The number of seed points (same as number of cells) in the Voronoi diagram
# - The maximum number of allowed iterations in Lloyd's algorithm
# - A tolerance on the maximum distance between a cell centroid and seed point.
diag = ot.LaguerreDiagram.BuildCentroidal(bbox,numPts,100,0.001)

# Plot the resulting centroidal Voronoi diagram
fig, ax = plt.subplots()
ot.PlotDiagram(diag, ax)

plt.savefig('CentroidalVoronoiDiagram.pdf')
plt.show()
