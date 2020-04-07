import pysdot as ot
import numpy as np
import matplotlib.pyplot as plt

numPts = 100

xbnds = [0.0,1.0] # minimum and maximum x values
ybnds = [0.0,1.0] # minimum and maximum y values
Ns = [50,50]

bbox = ot.BoundingBox(xbnds[0],xbnds[1],ybnds[0],ybnds[1])

grid = ot.RegularGrid(bbox, Ns[0], Ns[1])
dens = np.ones(Ns)
for i in range(Ns[0]):
    for j in range(Ns[1]):
        pt = grid.Center(i,j)
        dens[i,j] = np.exp(-30.0*( (pt[0]-0.5)**2 + (pt[1]-0.5)**2))

dist = ot.DiscretizedDistribution(grid,dens)

# Construct the Centroidal Voronoi diagram.  This function uses Lloyd's algorithm
# with latin hypercube samples as initial points (https://en.wikipedia.org/wiki/Lloyd%27s_algorithm)
# Arguments to BuildCentroidal are:
# - The bounding box
# - The number of seed points (same as number of cells) in the Voronoi diagram
# - The maximum number of allowed iterations in Lloyd's algorithm
# - A tolerance on the maximum distance between a cell centroid and seed point.
diag = ot.LaguerreDiagram.BuildCentroidal(bbox,numPts,1000,0.001,dist)

areas = diag.Areas(dist)

# Plot the resulting centroidal Voronoi diagram
fig, axs = plt.subplots(ncols=2,figsize=(14,6))
ot.PlotDiagram(diag, axs[0], distribution=dist, cell_colors=areas)
axs[0].set_title('Weighted CVD')

axs[1].imshow(dens.T,extent=[xbnds[0],xbnds[1],ybnds[0],ybnds[1]],origin='lower',alpha=0.8)
axs[1].set_title('Density')

plt.show()
