import pysdot as ot
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import powerlaw

numPts = 600

xbnds = [0.0,1.0] # minimum and maximum x values
ybnds = [0.0,1.0] # minimum and maximum y values
Ns = [1,1]

bbox = ot.BoundingBox(xbnds[0],xbnds[1],ybnds[0],ybnds[1])

grid = ot.RegularGrid(bbox, Ns[0], Ns[1])
dens = np.ones(Ns)
dist = ot.DiscretizedDistribution(grid,dens)


# The nugget is meant to prevent probabilities from getting too small and causing issues in the SDOT calculations
powExp = 0.3
nugget = 0.1
rv = powerlaw(powExp)
ptProbs = nugget+rv.rvs(size=numPts)

# Construct the Centroidal Power diagram.  This function uses Lloyd's algorithm
# with latin hypercube samples as initial points (https://en.wikipedia.org/wiki/Lloyd%27s_algorithm)
# Arguments to BuildCentroidal are:
# - The discretized distribution dist
# - The probabilities of the points
# - The maximum number of allowed iterations in Lloyd's algorithm
# - A tolerance on the maximum distance between a cell centroid and seed point.
diag = ot.SemidiscreteOT.BuildCentroidal(dist,ptProbs,10,0.001)

areas = diag.Areas(dist)

# Plot the resulting centroidal Voronoi diagram
fig, axs = plt.subplots(ncols=2,figsize=(14,6))
ot.PlotDiagram(diag, axs[0], distribution=dist, cell_colors=areas)
axs[0].set_title('Centroidal Power Diagram')

axs[1].hist(ptProbs,density=True)
v = np.linspace(0,1,150)
axs[1].plot(v+nugget, rv.pdf(v), linewidth=3)

axs[1].set_xlabel('Volume')
axs[1].set_title('Particle Volume Distribution')
axs[1].set_ylabel('Density')

plt.show()
